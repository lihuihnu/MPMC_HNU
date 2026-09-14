# Standard gRPC status schema

`status.proto` is an unmodified copy from Google APIs, commit
`ebd1d23ac613b177828dad42ad8dfb13ba498279` (the latest change to this file when acquired):
https://github.com/googleapis/googleapis/blob/ebd1d23ac613b177828dad42ad8dfb13ba498279/google/rpc/status.proto

Copyright Google LLC; Apache License 2.0, included in `LICENSE`. No upstream
NOTICE applies to this file. This does not select a license for MPMC_HNU itself.

Used solely for interoperable gRPC rich-error framing. Its `google.protobuf.Any`
import comes from the project's existing pinned Protobuf dependency. C++ and ES
bindings are generated during the build, never maintained by hand.
