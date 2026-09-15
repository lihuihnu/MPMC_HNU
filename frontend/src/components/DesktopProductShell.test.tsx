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
  it('keeps PT Flash as the default and does not render the Expert workspace', () => {
    const markup = html('pt');
    expect(markup).toContain('data-desktop-product-shell="true"');
    expect(markup).toContain('data-product-mode="pt"');
    expect(markup).toContain('aria-pressed="true">PT Flash');
    expect(markup).toContain('Model-neutral PT Flash');
    expect(markup).not.toContain('data-expert-workbench');
    expect(markup).not.toContain('PR76 model workspace');
  });

  it('opens the local editable workbench directly without account/session UI', () => {
    const markup = html('expert');
    expect(markup).toContain('data-expert-workbench="ready"');
    expect(markup).toContain('PR76 model workspace');
    expect(markup).toContain('Create custom PR76 model');
    expect(markup).toContain('No Expert model created yet');
    expect(markup).toContain('no account or login is required');
    expect(markup).not.toContain('authenticated model session');
    expect(markup).not.toContain('Reconnect Expert session');
  });
});
