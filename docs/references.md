# 项目参考资料

[返回项目入口](../README.md)。本页集中保留项目级工程参考与原始文献入口；具体采用的公式、参数、版本、许可和核验范围以各模块专题文档为准。列出参考不代表引入依赖或完成全文验证。工程软件的版本与能力须在实际引入或升级时重新核查。

- [E1] [CMake Presets 官方文档](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html)。
- [E2] [React：从头构建应用](https://react.dev/learn/build-a-react-app-from-scratch)；[Vite 官方指南](https://vite.dev/guide/)。
- [E3] [gRPC 核心概念](https://grpc.io/docs/what-is-grpc/core-concepts/)；[Protocol Buffers 语言指南](https://protobuf.dev/programming-guides/proto3/)。
- [E4] [gRPC-Web 官方实现与流支持说明](https://github.com/grpc/grpc-web#streaming-support)。
- [E5] [vtk.js 官方文档](https://kitware.github.io/vtk-js/docs/)。
- [E6] [CppAD 官方文档](https://coin-or.github.io/CppAD/)，用于 AD 模式、接口与独立核验方法参考，不预先绑定实现路线。
- [E7] [GitHub-hosted runners 官方文档](https://docs.github.com/actions/reference/runners/github-hosted-runners)。
- [E8] [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)；[OPM 模块说明](https://opm-project.org/?page_id=274)；[DuMux 官方文档](https://dumux.org/docs/doxygen/master/)。借鉴规范、分层与科学软件工程，复用代码前另行审计许可。
- [E9] [Connect for Web 官方生成与客户端文档](https://connectrpc.com/docs/web/generating-code/)；当前浏览器适配使用其 `createGrpcWebTransport` binary gRPC-Web 路径与 Protobuf-ES 生成描述符。
- [E10] [gRPC C++ `ServerBuilder`](https://grpc.github.io/grpc/cpp/classgrpc_1_1_server_builder.html)、[`ServerContext`](https://grpc.github.io/grpc/cpp/classgrpc_1_1_server_context.html)、[`ResourceQuota`](https://grpc.github.io/grpc/cpp/classgrpc_1_1_resource_quota.html)、[`SslServerCredentialsOptions`](https://grpc.github.io/grpc/cpp/structgrpc_1_1_ssl_server_credentials_options.html) 与 [authentication guide](https://grpc.io/docs/guides/auth/) 官方接口，以及 Envoy 官方 [`grpc_web`](https://www.envoyproxy.io/docs/envoy/latest/configuration/http/http_filters/grpc_web_filter)、[`cors`](https://www.envoyproxy.io/docs/envoy/latest/configuration/http/http_filters/cors_filter)、[`buffer`](https://www.envoyproxy.io/docs/envoy/latest/configuration/http/http_filters/buffer_filter)、[TLS transport](https://www.envoyproxy.io/docs/envoy/latest/api-v3/extensions/transport_sockets/tls/v3/tls.proto) 与 [gRPC health checking](https://www.envoyproxy.io/docs/envoy/latest/intro/arch_overview/upstream/health_checking) 文档；用于 process edge 的 message、memory、deadline/cancel、mTLS、health 与 browser transport policy。
- [S1] Peng, D.-Y.; Robinson, D. B. (1976). *A New Two-Constant Equation of State*. DOI: [10.1021/i160057a011](https://doi.org/10.1021/i160057a011)。
- [S2] Søreide, I.; Whitson, C. H. (1992). *Peng-Robinson predictions for hydrocarbons, CO2, N2, and H2S with pure water and NaCl brine*. DOI: [10.1016/0378-3812(92)85105-H](https://doi.org/10.1016/0378-3812(92)85105-H)。对应当前已实现的 SW92 corrected-original thermodynamics profile；具体已核验公式、勘误、参数语义与实现边界见 [SW92 thermodynamics 文档](../modules/thermodynamics/sw92.md)。
- [S3] Kontogeorgis, G. M.; Voutsas, E. C.; Yakoumis, I. V.; Tassios, D. P. (1996). *An Equation of State for Associating Fluids*. DOI: [10.1021/ie9600203](https://doi.org/10.1021/ie9600203)。
- [S4] Michelsen, M. L. (1982). *The isothermal flash problem. Part I. Stability*. DOI: [10.1016/0378-3812(82)85001-2](https://doi.org/10.1016/0378-3812(82)85001-2)。
- [S5] Michelsen, M. L. (1982). *The isothermal flash problem. Part II. Phase-split calculation*. DOI: [10.1016/0378-3812(82)85002-4](https://doi.org/10.1016/0378-3812(82)85002-4)。
- [S6] Xu, G.; Haynes, W. D.; Stadtherr, M. A. (2005). *Reliable Phase Stability Analysis for Asymmetric Models*. DOI: [10.1016/j.fluid.2005.06.016](https://doi.org/10.1016/j.fluid.2005.06.016)。当前 Xu-style 路径采用其 asymmetric lower-envelope/common-tangent formulation 构建 Gate 3A 与 maximum-two-phase Gate 3B；实现使用有限 multistart，不继承论文 interval-analysis 的全局保证。
