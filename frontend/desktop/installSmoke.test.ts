import { describe, expect, it } from 'vitest';

import {
  installedSmokeRequested,
  parseInstalledSmokePlan,
  PT_DESKTOP_INSTALL_SMOKE_PLAN_CONVENTION,
} from './installSmoke';

describe('installed desktop smoke contract', () => {
  it('is opt-in only through the dedicated command-line switch', () => {
    expect(installedSmokeRequested(['electron', 'app'])).toBe(false);
    expect(
      installedSmokeRequested(['electron', 'app', '--mpmc-install-smoke']),
    ).toBe(true);
  });

  it('parses only the versioned model-neutral request envelope', () => {
    const plan = parseInstalledSmokePlan({
      convention: PT_DESKTOP_INSTALL_SMOKE_PLAN_CONVENTION,
      expectedBackendCount: 1,
      requests: [
        {
          configuredBackendId: 'backend.fixture',
          pressurePa: 1e6,
          temperatureK: 300,
          feed: [{ componentId: 'component.fixture', moleFraction: 1 }],
        },
      ],
    });
    expect(plan.expectedBackendCount).toBe(1);
    expect(plan.requests[0]?.configuredBackendId).toBe('backend.fixture');
    expect(() =>
      parseInstalledSmokePlan({
        convention: 'wrong',
        expectedBackendCount: 1,
        requests: [],
      }),
    ).toThrow(/convention/u);
  });
});
