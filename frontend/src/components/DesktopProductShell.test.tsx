import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it, vi } from 'vitest';

import type { ExpertModelOwner, ExpertOwnedModel } from '../api/expertModelOwner';
import { unconfiguredFlashClient } from '../api/flashClient';
import { DesktopProductShellView } from './DesktopProductShell';

function owned(): ExpertOwnedModel {
  throw new Error('The static product-shell test must not create a model.');
}
function owner(): ExpertModelOwner {
  return { create: vi.fn(async () => owned()) };
}
function html(mode: 'pt' | 'expert') {
  return renderToStaticMarkup(
    <DesktopProductShellView
      flashClient={unconfiguredFlashClient}
      expertOwner={owner()}
      mode={mode}
      onPtMode={() => {}}
      onExpertMode={() => {}}
    />,
  );
}

describe('Electron desktop product shell', () => {
  it('presents the no-login classic PR workspace as the primary product mode', () => {
    const markup = html('expert');
    expect(markup).toContain('data-desktop-product-shell="true"');
    expect(markup).toContain('data-product-mode="expert"');
    expect(markup).toContain('aria-pressed="true">Classic PR');
    expect(markup).toContain('data-expert-workbench="ready"');
    expect(markup).toContain('Peng–Robinson flash');
    expect(markup).toContain('Create custom PR76 model');
    expect(markup).toContain('No Expert model created yet');
    expect(markup).toContain('no account or login required');
    expect(markup).not.toContain('authenticated model session');
    expect(markup).not.toContain('Reconnect Expert session');
  });

  it('keeps the pre-existing configured PT surface as a secondary compatibility mode', () => {
    const markup = html('pt');
    expect(markup).toContain('data-product-mode="pt"');
    expect(markup).toContain('aria-pressed="true">Configured PT');
    expect(markup).toContain('Model-neutral PT Flash');
    expect(markup).not.toContain('data-expert-workbench');
  });
});
