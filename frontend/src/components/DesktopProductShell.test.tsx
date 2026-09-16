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
function html(
  mode: 'pt' | 'expert',
  labels: { environmentLabel?: string; workspaceAccessLabel?: string } = {},
) {
  return renderToStaticMarkup(
    <DesktopProductShellView
      flashClient={unconfiguredFlashClient}
      expertOwner={owner()}
      mode={mode}
      onPtMode={() => {}}
      onExpertMode={() => {}}
      {...labels}
    />,
  );
}

describe('shared Classic PR product shell', () => {
  it('preserves the existing local no-login labels by default', () => {
    const markup = html('expert');
    expect(markup).toContain('data-desktop-product-shell="true"');
    expect(markup).toContain('data-product-mode="expert"');
    expect(markup).toContain('aria-pressed="true">Classic PR');
    expect(markup).toContain('data-expert-workbench="ready"');
    expect(markup).toContain('Peng–Robinson flash');
    expect(markup).toContain('Define PR fluid');
    expect(markup).toContain('Create a PR fluid model');
    expect(markup).toContain('no account or login required');
    expect(markup).toContain('Local calculation · no login');
    expect(markup).not.toContain('authenticated model session');
    expect(markup).not.toContain('Reconnect Expert session');
  });

  it('supports hosted authenticated labels without changing the shared Expert workspace', () => {
    const markup = html('expert', {
      environmentLabel: 'Hosted calculation · authenticated model session',
      workspaceAccessLabel: 'Hosted calculation · authenticated session',
    });
    expect(markup).toContain('Hosted calculation · authenticated model session');
    expect(markup).toContain('Hosted calculation · authenticated session');
    expect(markup).toContain('data-expert-workspace="pr76"');
    expect(markup).toContain('Peng–Robinson flash');
    expect(markup).not.toContain('Local calculation · no login');
    expect(markup).not.toContain('no account or login required');
  });

  it('keeps the pre-existing configured PT surface as a secondary compatibility mode', () => {
    const markup = html('pt');
    expect(markup).toContain('data-product-mode="pt"');
    expect(markup).toContain('aria-pressed="true">Configured PT');
    expect(markup).toContain('Model-neutral PT Flash');
    expect(markup).not.toContain('data-expert-workbench');
  });
});
