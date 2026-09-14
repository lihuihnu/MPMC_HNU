package org.mpmc.ptandroid;

import android.app.Activity;
import android.os.Bundle;
import android.util.Log;

public final class SmokeActivity extends Activity {
    private static final String TAG = "MPMC_PT_ANDROID";

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        Thread worker = new Thread(
            () -> {
                String result;
                try {
                    result = NativeBridge.runSmoke();
                } catch (Throwable error) {
                    result = "ANDROID_JNI_PT_SERVICE_FAIL java="
                        + error.getClass().getSimpleName();
                }

                if (result.startsWith("ANDROID_JNI_PT_SERVICE_OK")) {
                    Log.i(TAG, result);
                } else {
                    Log.e(TAG, result);
                }
                runOnUiThread(this::finishAndRemoveTask);
            },
            "mpmc-pt-native-smoke"
        );
        worker.start();
    }
}
