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
}

export function DesktopProductShellView({
  flashClient,
  expertOwner,
  mode,
  onPtMode,
  onExpertMode,
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
          <span>Local desktop calculation · no account or login required</span>
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
          <ExpertPr76Workspace owner={expertOwner} />
        </div>
      ) : null}
      {mode === 'pt' ? <App client={flashClient} /> : null}
    </>
  );
}

export interface DesktopProductShellProps {
  flashClient: FlashClient;
  expertOwner: ExpertModelOwner;
}

/**
 * Electron-only shell. Classic PR is the product entry point and requires no
 * account/session action; the pre-existing configured-PT surface remains a
 * secondary compatibility path.
 */
export function DesktopProductShell({ flashClient, expertOwner }: DesktopProductShellProps) {
  const [mode, setMode] = useState<DesktopProductMode>('expert');
  return (
    <DesktopProductShellView
      flashClient={flashClient}
      expertOwner={expertOwner}
      mode={mode}
      onPtMode={() => setMode('pt')}
      onExpertMode={() => setMode('expert')}
    />
  );
}
