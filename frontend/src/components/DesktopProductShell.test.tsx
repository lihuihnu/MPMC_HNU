import { renderToStaticMarkup } from 'react-dom/server';
import { Code } from '@connectrpc/connect';
import { describe, expect, it, vi } from 'vitest';

import type { ExpertModelOwner, ExpertOwnedModel } from '../api/expertModelOwner';
import { unconfiguredFlashClient } from '../api/flashClient';
import { RendererModelError } from '../api/rendererModelClient';
import {
  DesktopProductShellView,
  desktopExpertFailure,
  type DesktopExpertSessionState,
} from './DesktopProductShell';

function owned(): ExpertOwnedModel {
  throw new Error('The static product-shell test must not create a model.');
}
function owner(): ExpertModelOwner {
  return { create: vi.fn(async () => owned()) };
}
function html(mode: 'pt' | 'expert', expertSession: DesktopExpertSessionState) {
  return renderToStaticMarkup(
    <DesktopProductShellView
      flashClient={unconfiguredFlashClient}
      expertOwner={owner()}
      mode={mode}
      expertSession={expertSession}
      onPtMode={() => {}}
      onExpertMode={() => {}}
      onReconnect={() => {}}
    />,
  );
}

describe('Electron desktop product shell', () => {
  it('keeps PT Flash as the default and does not render/open Expert workspace', () => {
    const markup = html('pt', { status: 'idle' });
    expect(markup).toContain('data-desktop-product-shell="true"');
    expect(markup).toContain('data-product-mode="pt"');
    expect(markup).toContain('aria-pressed="true">PT Flash');
    expect(markup).toContain('Model-neutral PT Flash');
    expect(markup).not.toContain('data-expert-session');
    expect(markup).not.toContain('PR76 model workspace');
  });

  it('renders explicit loading and failure states without inventing a model', () => {
    const loading = html('expert', { status: 'loading' });
    expect(loading).toContain('data-expert-session="loading"');
    expect(loading).toContain('Opening authenticated model session');
    expect(loading).not.toContain('PR76 model workspace');

    const failed = html('expert', {
      status: 'failed',
      failure: { reason: 'renderer.reconnect_required', code: Code.FailedPrecondition },
    });
    expect(failed).toContain('data-expert-session="failed"');
    expect(failed).toContain('renderer.reconnect_required');
    expect(failed).toContain('Reconnect Expert session');
    expect(failed).toContain('Reconnection is never automatic');
  });

  it('renders the verified workspace only after the authenticated session is ready', () => {
    const markup = html('expert', { status: 'ready' });
    expect(markup).toContain('data-expert-session="ready"');
    expect(markup).toContain('PR76 model workspace');
    expect(markup).toContain('Create custom PR76 model');
    expect(markup).toContain('No Expert model created yet');
  });

  it('preserves typed session failures and sanitizes arbitrary exception text', () => {
    expect(desktopExpertFailure(
      new RendererModelError(Code.Unavailable, 'session.open_failed', 'ipc'),
    )).toEqual({ reason: 'session.open_failed', code: Code.Unavailable });
    const unknown = desktopExpertFailure(new Error('private-host-token'));
    expect(unknown).toEqual({ reason: 'desktop.expert_session_failed' });
    expect(JSON.stringify(unknown)).not.toContain('private-host-token');
  });
});
