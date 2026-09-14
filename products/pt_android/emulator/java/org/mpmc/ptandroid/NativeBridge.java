package org.mpmc.ptandroid;

public final class NativeBridge {
    static {
        System.loadLibrary("mpmc_pt_android_core");
    }

    private NativeBridge() {}

    public static native String runSmoke();
}
