package org.mpmc.ptandroid;

final class NativeBridge {
    static {
        System.loadLibrary("mpmc_pt_android_core");
    }

    private NativeBridge() {}

    static native String discoverJson();

    static native String solveJson(
            String configuredBackendId,
            double pressurePa,
            double temperatureK,
            String[] componentIds,
            double[] moleFractions);

    static native String[] modelApplyRecords(String requestId, String[] records);

    static native String[] modelSolveRecords(String requestId, String[] records);

    static native String[] modelReleaseRecords(String requestId);

    static native void modelCancel(String requestId);
}
