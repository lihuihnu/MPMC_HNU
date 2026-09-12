import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it } from 'vitest';

import { mapDiscoveryResponse } from '../api/ptWire';
import { wireDiscovery } from '../test/ptWireFixtures';
import { FlashForm } from './FlashForm';

describe('discovery-driven PT form', () => {
  it('renders backend and component identities only from discovery', () => {
    const backends = mapDiscoveryResponse(wireDiscovery()).backends;
    const html = renderToStaticMarkup(
      <FlashForm
        backends={backends}
        discoveryLoading={false}
        busy={false}
        onBackendChange={() => undefined}
        onCancel={() => undefined}
        onSubmit={() => undefined}
      />,
    );
    expect(html).toContain('pr76-default');
    expect(html).toContain('methane');
    expect(html).toContain('carbon-dioxide');
    expect(html).toContain('IDs and order come from discovery');
    expect(html).not.toContain('NaCl molality');
    expect(html).not.toContain('Add component');
  });
});
