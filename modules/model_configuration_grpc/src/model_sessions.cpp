#include <mpmc/model_configuration_grpc/model_sessions.hpp>
#include <google/protobuf/unknown_field_set.h>

#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <utility>

namespace mpmc::model_configuration_grpc {
namespace {
namespace wire = ::mpmc::model_configuration::v1;
namespace mc = ::mpmc::model_configuration;
using SC = grpc::StatusCode;
struct Failure { SC status; const char* code; };
[[noreturn]] void fail(SC status, const char* code) { throw Failure{status, code}; }
grpc::Status failure_status() noexcept {
    try {
        SC status = SC::INTERNAL;
        const char* code = "session.internal_failure";
        try { throw; }
        catch (const Failure& e) { status = e.status; code = e.code; }
        catch (const std::bad_alloc&) { status = SC::RESOURCE_EXHAUSTED; code = "session.memory_exhausted"; }
        catch (const mc::ModelRegistryError& e) {
            if (e.code() == mc::ModelRegistryErrorCode::entropy_unavailable) {
                status = SC::UNAVAILABLE; code = "session.entropy_unavailable";
            }
        }
        catch (...) {}
        wire::ModelServiceError detail;
        detail.set_wire_contract(std::string(model_wire_contract)); detail.set_code(code);
        return model_error_status(status, "model session request failed", detail);
    } catch (...) { return {SC::INTERNAL, "model session error mapping failed"}; }
}
struct Count {
    std::atomic<std::size_t>& count;
    Count(std::atomic<std::size_t>& value, std::size_t maximum) : count(value) {
        auto n = count.load(std::memory_order_relaxed);
        while (n < maximum) {
            if (count.compare_exchange_weak(n, n + 1, std::memory_order_acquire)) { return; }
        }
        fail(SC::RESOURCE_EXHAUSTED, "session.capacity_exceeded");
    }
    ~Count() { count.fetch_sub(1, std::memory_order_release); }
    Count(const Count&) = delete;
    Count& operator=(const Count&) = delete;
};
bool valid_session_id(std::string_view id) {
    if (id.size() != 84 || !id.starts_with("ms1_")) { return false; }
    return std::all_of(id.begin() + 4, id.end(), [](char c) {
        return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
    });
}
} // namespace
struct ModelSessionService::Impl {
    struct Entry {
        // Slot dies last, after adapter/registry and all retained handlers.
        Count slot;
        const std::string principal;
        mc::Pr76ModelRegistry registry;
        ModelGrpcServiceAdapter adapter;
        Entry(Impl& host, std::string owner)
            : slot(host.residents, host.policy.max_sessions), principal(std::move(owner)),
              registry(mc::Pr76ModelRegistryLimits{host.policy.max_models_per_session, {}, {}},
                       mc::ModelDataPolicy::ordinary, host.entropy),
              // Private adapter is reached only after this router authenticates
              // identity and session. No unauthenticated listener registers it.
              adapter(registry, [](const grpc::ServerContext&) { return true; }, host.policy.rpc) {}
    };
    ModelSessionIdentity identity;
    const ModelSessionLimits policy;
    mc::ModelHandleEntropySource entropy;
    std::atomic<std::size_t> residents{0}, requests{0}, solves{0};
    mutable std::mutex mutex;
    std::condition_variable changed;
    std::map<std::string, std::shared_ptr<Entry>, std::less<>> entries;
    std::uint64_t sequence{0};
    bool closed{false};
    Impl(ModelSessionIdentity auth, ModelSessionLimits limits, mc::ModelHandleEntropySource source)
        : identity(std::move(auth)), policy(limits), entropy(std::move(source)) {
        if (!identity || !entropy || !policy.rpc.structurally_valid() || policy.max_sessions == 0 ||
            policy.max_models_per_session == 0 || policy.max_session_lifetime.count() <= 0 ||
            policy.max_session_lifetime > std::chrono::hours(24)) {
            throw std::invalid_argument("invalid host model session policy");
        }
    }
    std::string authenticate(const grpc::ServerContext& context) const {
        auto principal = identity(context);
        if (principal.empty() || principal.size() > 64U * 1024U) {
            fail(SC::UNAUTHENTICATED, "session.unauthenticated");
        }
        return principal;
    }
    std::pair<std::string, std::shared_ptr<Entry>> create(std::string principal) {
        auto entry = std::make_shared<Entry>(*this, std::move(principal));
        const auto random = entropy(); // Capacity reserved before entropy/allocation.
        static constexpr char hex[] = "0123456789abcdef";
        std::string id("ms1_"); id.reserve(84);
        for (auto byte : random) { id += hex[byte >> 4U]; id += hex[byte & 15U]; }
        std::lock_guard lock(mutex);
        if (closed) { fail(SC::FAILED_PRECONDITION, "session.host_closed"); }
        if (sequence == std::numeric_limits<std::uint64_t>::max()) {
            fail(SC::RESOURCE_EXHAUSTED, "session.capacity_exceeded");
        }
        const auto serial = ++sequence;
        for (int shift = 60; shift >= 0; shift -= 4) { id += hex[(serial >> shift) & 15U]; }
        entries.emplace(id, entry);
        return {std::move(id), std::move(entry)};
    }
    void remove(const std::string& id, const std::shared_ptr<Entry>& entry) {
        {
            std::lock_guard lock(mutex);
            entries.erase(id); // Future lookup fails before registry closure.
        }
        entry->registry.close();
        changed.notify_all();
    }
    std::shared_ptr<Entry> lookup(const grpc::ServerContext& context) {
        const auto principal = authenticate(context);
        const auto range = context.client_metadata().equal_range(model_session_metadata);
        if (range.first == range.second || std::next(range.first) != range.second) {
            fail(SC::INVALID_ARGUMENT, "session.metadata_required");
        }
        const auto& raw = range.first->second;
        const std::string_view id(raw.data(), raw.size());
        if (!valid_session_id(id)) { fail(SC::INVALID_ARGUMENT, "session.invalid_id"); }
        std::lock_guard lock(mutex);
        if (closed) { fail(SC::FAILED_PRECONDITION, "session.host_closed"); }
        const auto found = entries.find(id);
        if (found == entries.end()) { fail(SC::NOT_FOUND, "session.not_found"); }
        if (found->second->principal != principal) { fail(SC::PERMISSION_DENIED, "session.wrong_identity"); }
        return found->second;
    }
    template <class Request, class Response, class Method>
    grpc::Status forward(grpc::ServerContext* context, const Request* request, Response* response,
                         Method method, bool solving = false) {
        response->Clear();
        try {
            const Count admitted(requests, policy.rpc.max_concurrent_requests);
            const auto entry = lookup(*context);
            if (solving) {
                const Count solve_admission(solves, policy.rpc.max_concurrent_solves);
                return (entry->adapter.*method)(context, request, response);
            }
            return (entry->adapter.*method)(context, request, response);
        } catch (...) { response->Clear(); return failure_status(); }
    }
};
ModelSessionService::ModelSessionService(ModelSessionIdentity identity, ModelSessionLimits limits,
    mc::ModelHandleEntropySource entropy)
    : impl_(std::make_unique<Impl>(std::move(identity), limits, std::move(entropy))) {}
ModelSessionService::~ModelSessionService() { close(); }
const ModelSessionLimits& ModelSessionService::limits() const noexcept { return impl_->policy; }
ModelSessionStatus ModelSessionService::status() const {
    std::lock_guard lock(impl_->mutex);
    return {impl_->entries.size(), impl_->residents.load(), impl_->requests.load(), impl_->closed};
}
void ModelSessionService::close() {
    decltype(impl_->entries) removed;
    {
        std::lock_guard lock(impl_->mutex);
        impl_->closed = true;
        removed.swap(impl_->entries);
    }
    for (auto& [id, entry] : removed) { (void)id; entry->registry.close(); }
    impl_->changed.notify_all();
}
void ModelSessionService::register_services(grpc::ServerBuilder& builder) {
    builder.RegisterService(static_cast<wire::ModelConfigurationService::Service*>(this));
    builder.RegisterService(static_cast<wire::ModelSessionService::Service*>(this));
}
grpc::Status ModelSessionService::OpenModelSession(grpc::ServerContext* context,
    const wire::OpenModelSessionRequest* request, grpc::ServerWriter<wire::ModelSessionOpened>* writer) {
    std::string id;
    std::shared_ptr<Impl::Entry> entry;
    try {
        auto principal = impl_->authenticate(*context);
        const auto now = std::chrono::system_clock::now();
        if (context->deadline() == std::chrono::system_clock::time_point::max()) {
            fail(SC::INVALID_ARGUMENT, "session.deadline_required");
        }
        if (context->deadline() <= now || context->IsCancelled()) {
            fail(SC::CANCELLED, "session.cancelled");
        }
        if (context->deadline() - now > impl_->policy.max_session_lifetime) {
            fail(SC::INVALID_ARGUMENT, "session.deadline_required");
        }
        if (request->ByteSizeLong() > impl_->policy.rpc.max_request_bytes ||
            request->wire_contract() != model_session_contract ||
            request->GetReflection()->GetUnknownFields(*request).field_count() != 0) {
            fail(SC::INVALID_ARGUMENT, "session.invalid_request");
        }
        auto created = impl_->create(std::move(principal));
        id = std::move(created.first); entry = std::move(created.second);
        wire::ModelSessionOpened opened;
        opened.set_wire_contract(std::string(model_session_contract)); opened.set_session_id(id);
        if (opened.ByteSizeLong() > impl_->policy.rpc.max_response_bytes || !writer->Write(opened)) {
            fail(SC::CANCELLED, "session.publication_failed");
        }
        // gRPC marks cancellation on stream failure/deadline. A finite deadline
        // also bounds retention when a dead network is not immediately detected.
        std::unique_lock lock(impl_->mutex);
        while (!impl_->closed && !context->IsCancelled() &&
               std::chrono::system_clock::now() < context->deadline()) {
            impl_->changed.wait_for(lock, std::chrono::milliseconds(20));
        }
        lock.unlock();
        impl_->remove(id, entry);
        return {SC::CANCELLED, "model session closed"};
    } catch (...) {
        const auto status = failure_status();
        if (entry) { impl_->remove(id, entry); }
        return status;
    }
}
grpc::Status ModelSessionService::CreateModel(grpc::ServerContext* c, const wire::CreateModelRequest* r, wire::CreateModelResponse* o) {
    return impl_->forward(c, r, o, &ModelGrpcServiceAdapter::CreateModel);
}
grpc::Status ModelSessionService::DescribeModel(grpc::ServerContext* c, const wire::ModelHandleRequest* r, wire::DescribeModelResponse* o) {
    return impl_->forward(c, r, o, &ModelGrpcServiceAdapter::DescribeModel);
}
grpc::Status ModelSessionService::SolveModel(grpc::ServerContext* c, const wire::SolveModelRequest* r, wire::SolveModelResponse* o) {
    return impl_->forward(c, r, o, &ModelGrpcServiceAdapter::SolveModel, true);
}
grpc::Status ModelSessionService::ReleaseModel(grpc::ServerContext* c, const wire::ModelHandleRequest* r, wire::ReleaseModelResponse* o) {
    return impl_->forward(c, r, o, &ModelGrpcServiceAdapter::ReleaseModel);
}
} // namespace mpmc::model_configuration_grpc
