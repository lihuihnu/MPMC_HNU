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
}
