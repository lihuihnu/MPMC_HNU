package org.mpmc.ptandroid;

import com.getcapacitor.JSObject;
import com.getcapacitor.Plugin;
import com.getcapacitor.PluginCall;
import com.getcapacitor.PluginMethod;
import com.getcapacitor.annotation.CapacitorPlugin;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

@CapacitorPlugin(name = "MpmcPt")
public final class MpmcPtPlugin extends Plugin {
    private static final ExecutorService NATIVE_EXECUTOR =
            Executors.newSingleThreadExecutor(runnable -> {
                Thread thread = new Thread(runnable, "mpmc-pt-native");
                thread.setDaemon(true);
                return thread;
            });

    @PluginMethod
    public void discover(PluginCall call) {
        NATIVE_EXECUTOR.execute(() -> resolveJson(call, NativeBridge.discoverJson()));
    }

    @PluginMethod
    public void solve(PluginCall call) {
        final String requestJson = call.getString("requestJson");
        if (requestJson == null || requestJson.isEmpty()) {
            call.reject("requestJson is required.", "MPMC_PT_ANDROID_INVALID_REQUEST");
            return;
        }

        NATIVE_EXECUTOR.execute(() -> {
            try {
                JSONObject request = new JSONObject(requestJson);
                String configuredBackendId = request.getString("configuredBackendId");
                double pressurePa = request.getDouble("pressurePa");
                double temperatureK = request.getDouble("temperatureK");
                JSONArray feed = request.getJSONArray("feed");
                String[] componentIds = new String[feed.length()];
                double[] moleFractions = new double[feed.length()];
                for (int index = 0; index < feed.length(); ++index) {
                    JSONObject component = feed.getJSONObject(index);
                    componentIds[index] = component.getString("componentId");
                    moleFractions[index] = component.getDouble("moleFraction");
                }
                resolveJson(
                        call,
                        NativeBridge.solveJson(
                                configuredBackendId,
                                pressurePa,
                                temperatureK,
                                componentIds,
                                moleFractions));
            } catch (JSONException error) {
                rejectOnUiThread(
                        call,
                        "Android PT request JSON is malformed: " + error.getMessage(),
                        "MPMC_PT_ANDROID_INVALID_REQUEST");
            } catch (RuntimeException error) {
                rejectOnUiThread(
                        call,
                        "Android PT native bridge failed: " + error.getMessage(),
                        "MPMC_PT_ANDROID_BRIDGE_FAILURE");
            }
        });
    }

    private void resolveJson(PluginCall call, String json) {
        if (json == null || json.isEmpty()) {
            rejectOnUiThread(
                    call,
                    "Android PT native bridge returned an empty payload.",
                    "MPMC_PT_ANDROID_EMPTY_REPLY");
            return;
        }
        JSObject result = new JSObject();
        result.put("json", json);
        getActivity().runOnUiThread(() -> call.resolve(result));
    }

    private void rejectOnUiThread(PluginCall call, String message, String code) {
        getActivity().runOnUiThread(() -> call.reject(message, code));
    }
}
