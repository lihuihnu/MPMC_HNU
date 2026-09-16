import { useState } from 'react';

import { App } from '../App';
import type { ExpertModelOwner } from '../api/expertModelOwner';
import type { FlashClient } from '../api/flashClient';
import { ExpertPr76Workspace } from './ExpertPr76Workspace';

export type DesktopProductMode = 'pt' | 'expert';

export interface DesktopProductShellViewProps {
  flashClient: FlashClient;
  expertOwner: ExpertModelOwner;
  mode: DesktopProductMode;
  onPtMode(): void;
  onExpertMode(): void;
  environmentLabel?: string;
  workspaceAccessLabel?: string;
}

export function DesktopProductShellView({
  flashClient,
  expertOwner,
  mode,
  onPtMode,
  onExpertMode,
  environmentLabel = 'Local desktop calculation · no account or login required',
  workspaceAccessLabel = 'Local calculation · no login',
}: DesktopProductShellViewProps) {
  return (
    <>
      <nav
        className="desktop-product-mode"
        data-desktop-product-shell="true"
        aria-label="Flash workspace mode"
      >
        <div>
          <strong>Flash workspace</strong>
          <span>{environmentLabel}</span>
        </div>
        <div className="submit-actions">
          <button
            type="button"
            className={mode === 'expert' ? 'primary-button' : 'secondary-button cancel-button'}
            data-product-mode="expert"
            aria-pressed={mode === 'expert'}
            onClick={onExpertMode}
          >
            Classic PR
          </button>
          <button
            type="button"
            className={mode === 'pt' ? 'primary-button' : 'secondary-button cancel-button'}
            data-product-mode="pt"
            aria-pressed={mode === 'pt'}
            onClick={onPtMode}
          >
            Configured PT
          </button>
        </div>
      </nav>

      {mode === 'expert' ? (
        <div data-expert-workbench="ready">
          <ExpertPr76Workspace owner={expertOwner} accessLabel={workspaceAccessLabel} />
        </div>
      ) : null}
      {mode === 'pt' ? <App client={flashClient} /> : null}
    </>
  );
}

export interface DesktopProductShellProps {
  flashClient: FlashClient;
  expertOwner: ExpertModelOwner;
  environmentLabel?: string;
  workspaceAccessLabel?: string;
}

/**
 * Shared Classic PR product shell. Existing local products use the default
 * no-login labels; hosted callers may override presentation text without
 * changing ownership or scientific semantics.
 */
export function DesktopProductShell({
  flashClient,
  expertOwner,
  environmentLabel,
  workspaceAccessLabel,
}: DesktopProductShellProps) {
  const [mode, setMode] = useState<DesktopProductMode>('expert');
  return (
    <DesktopProductShellView
      flashClient={flashClient}
      expertOwner={expertOwner}
      mode={mode}
      onPtMode={() => setMode('pt')}
      onExpertMode={() => setMode('expert')}
      {...(environmentLabel === undefined ? {} : { environmentLabel })}
      {...(workspaceAccessLabel === undefined ? {} : { workspaceAccessLabel })}
    />
  );
}
