package org.mpmc.ptandroid;

import android.os.Bundle;
import com.getcapacitor.BridgeActivity;

public final class MainActivity extends BridgeActivity {
    @Override
    public void onCreate(Bundle savedInstanceState) {
        registerPlugin(MpmcPtPlugin.class);
        super.onCreate(savedInstanceState);
    }
}
