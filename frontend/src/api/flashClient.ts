import type { ProfileCPtRequest, ProfileCPtResponse } from '../domain/flash';

export interface FlashClient {
  readonly configured: boolean;
  solveProfileCPt(request: ProfileCPtRequest): Promise<ProfileCPtResponse>;
}

export class BackendNotConfiguredError extends Error {
  constructor() {
    super('The frontend build has no configured MPMC_HNU compute service.');
    this.name = 'BackendNotConfiguredError';
  }
}

export const unconfiguredFlashClient: FlashClient = {
  configured: false,
  async solveProfileCPt(_request: ProfileCPtRequest): Promise<ProfileCPtResponse> {
    throw new BackendNotConfiguredError();
  },
};
