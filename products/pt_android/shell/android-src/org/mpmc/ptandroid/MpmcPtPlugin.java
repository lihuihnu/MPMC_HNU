package org.mpmc.ptandroid;

import android.util.Log;
import com.getcapacitor.JSObject;
import com.getcapacitor.Plugin;
import com.getcapacitor.PluginCall;
import com.getcapacitor.PluginMethod;
import com.getcapacitor.annotation.CapacitorPlugin;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

@CapacitorPlugin(name = "MpmcPt")
public final class MpmcPtPlugin extends Plugin {
    private static final String LOG_TAG = "MPMC_PT_ANDROID";
    private static final int MAX_MODEL_REQUEST_BYTES = 64 * 1024;
    private static final ExecutorService NATIVE_EXECUTOR =
            Executors.newSingleThreadExecutor(runnable -> {
                Thread thread = new Thread(runnable, "mpmc-pt-native");
                thread.setDaemon(true);
                return thread;
            });

    @PluginMethod
    public void discover(PluginCall call) {
        NATIVE_EXECUTOR.execute(() -> {
            String json = NativeBridge.discoverJson();
            logDiscoveryEvidence(json);
            resolveJson(call, json);
        });
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
                String json = NativeBridge.solveJson(
                        configuredBackendId,
                        pressurePa,
                        temperatureK,
                        componentIds,
                        moleFractions);
                logSolveEvidence(configuredBackendId, json);
                resolveJson(call, json);
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

    @PluginMethod
    public void modelApply(PluginCall call) {
        invokeModel(call, "apply");
    }

    @PluginMethod
    public void modelSolve(PluginCall call) {
        invokeModel(call, "solve");
    }

    @PluginMethod
    public void modelRelease(PluginCall call) {
        final String requestId = requiredRequestId(call);
        if (requestId == null) return;
        NATIVE_EXECUTOR.execute(() -> {
            String[] records = NativeBridge.modelReleaseRecords(requestId);
            logModelEvidence("RELEASE", records);
            resolveRecords(call, records);
        });
    }

    @PluginMethod
    public void modelCancel(PluginCall call) {
        final String requestId = requiredRequestId(call);
        if (requestId == null) return;
        // Cancellation must not queue behind a long native solve/apply. Native
        // ownership state is mutex protected; the worker safely observes this mark.
        NativeBridge.modelCancel(requestId);
        call.resolve();
    }

    private void invokeModel(PluginCall call, String operation) {
        final String requestId = requiredRequestId(call);
        if (requestId == null) return;
        final String recordsJson = call.getString("recordsJson");
        if (recordsJson == null || recordsJson.isEmpty()
                || recordsJson.getBytes(StandardCharsets.UTF_8).length > MAX_MODEL_REQUEST_BYTES) {
            call.reject("recordsJson is required and bounded.", "MPMC_PT_ANDROID_INVALID_REQUEST");
            return;
        }
        final String[] records;
        try {
            JSONArray input = new JSONArray(recordsJson);
            records = new String[input.length()];
            for (int index = 0; index < input.length(); ++index) {
                records[index] = input.getString(index);
                if (records[index].indexOf('\0') >= 0) {
                    throw new JSONException("embedded NUL");
                }
            }
        } catch (JSONException error) {
            call.reject("Android model records are malformed.", "MPMC_PT_ANDROID_INVALID_REQUEST");
            return;
        }
        NATIVE_EXECUTOR.execute(() -> {
            String[] reply = "apply".equals(operation)
                    ? NativeBridge.modelApplyRecords(requestId, records)
                    : NativeBridge.modelSolveRecords(requestId, records);
            logModelEvidence("apply".equals(operation) ? "APPLY" : "SOLVE", reply);
            resolveRecords(call, reply);
        });
    }

    private String requiredRequestId(PluginCall call) {
        String requestId = call.getString("requestId");
        if (requestId == null || requestId.isEmpty() || requestId.length() > 128
                || requestId.indexOf('\0') >= 0) {
            call.reject("requestId is invalid.", "MPMC_PT_ANDROID_INVALID_REQUEST");
            return null;
        }
        return requestId;
    }

    // Product evidence never logs model/component parameters or request state.
    private static void logModelEvidence(String operation, String[] records) {
        if (records != null && records.length > 0 && records[0].startsWith("ok")) {
            Log.i(LOG_TAG, "ANDROID_CLASSIC_PR_" + operation + "_OK");
        }
    }

    // Product Shell v1 emits bounded operational evidence on all builds. These
    // markers contain only configured backend identity and phase count; request
    // state, compositions, thermodynamic parameters, and credentials are never
    // logged. Avoid depending on generated BuildConfig so the app-local plugin
    // remains template-neutral across Capacitor Android revisions.
    private static void logDiscoveryEvidence(String json) {
        try {
            JSONObject payload = new JSONObject(json);
            JSONArray backends = payload.optJSONArray("backends");
            if (backends != null) {
                Log.i(LOG_TAG, "ANDROID_PRODUCT_SHELL_DISCOVERY_OK backends=" + backends.length());
            }
        } catch (JSONException ignored) {
            Log.w(LOG_TAG, "ANDROID_PRODUCT_SHELL_DISCOVERY_REPLY_INVALID_JSON");
        }
    }

    private static void logSolveEvidence(String configuredBackendId, String json) {
        try {
            JSONObject payload = new JSONObject(json);
            JSONObject response = payload.optJSONObject("response");
            if (response == null || !"result".equals(response.optString("kind"))) return;
            JSONObject result = response.optJSONObject("result");
            JSONArray phases = result == null ? null : result.optJSONArray("phases");
            if (phases != null) {
                Log.i(
                        LOG_TAG,
                        "ANDROID_PRODUCT_SHELL_SOLVE_OK backend="
                                + configuredBackendId
                                + " phases="
                                + phases.length());
            }
        } catch (JSONException ignored) {
            Log.w(LOG_TAG, "ANDROID_PRODUCT_SHELL_SOLVE_REPLY_INVALID_JSON");
        }
    }

    private void resolveJson(PluginCall call, String json) {
        if (json == null || json.isEmpty()) {
            rejectOnUiThread(call, "Android PT native bridge returned an empty payload.", "MPMC_PT_ANDROID_EMPTY_REPLY");
            return;
        }
        JSObject result = new JSObject();
        result.put("json", json);
        getActivity().runOnUiThread(() -> call.resolve(result));
    }

    private void resolveRecords(PluginCall call, String[] records) {
        if (records == null || records.length == 0) {
            rejectOnUiThread(call, "Android model bridge returned an empty payload.", "MPMC_PT_ANDROID_EMPTY_REPLY");
            return;
        }
        JSONArray values = new JSONArray();
        for (String record : records) values.put(record);
        JSObject result = new JSObject();
        result.put("records", values);
        getActivity().runOnUiThread(() -> call.resolve(result));
    }

    private void rejectOnUiThread(PluginCall call, String message, String code) {
        getActivity().runOnUiThread(() -> call.reject(message, code));
    }
}
