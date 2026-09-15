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
      <nav className="desktop-product-mode" data-desktop-product-shell="true" aria-label="Desktop product mode">
        <div>
          <strong>Desktop workspace</strong>
          <span>PT Flash and the editable PR76 workbench run locally; no account or login is required.</span>
        </div>
        <div className="submit-actions">
          <button
            type="button"
            className={mode === 'pt' ? 'primary-button' : 'secondary-button cancel-button'}
            data-product-mode="pt"
            aria-pressed={mode === 'pt'}
            onClick={onPtMode}
          >
            PT Flash
          </button>
          <button
            type="button"
            className={mode === 'expert' ? 'primary-button' : 'secondary-button cancel-button'}
            data-product-mode="expert"
            aria-pressed={mode === 'expert'}
            onClick={onExpertMode}
          >
            PR76 Expert
          </button>
        </div>
      </nav>

      {mode === 'pt' ? <App client={flashClient} /> : null}
      {mode === 'expert' ? (
        <div data-expert-workbench="ready">
          <ExpertPr76Workspace owner={expertOwner} />
        </div>
      ) : null}
    </>
  );
}

export interface DesktopProductShellProps {
  flashClient: FlashClient;
  expertOwner: ExpertModelOwner;
}

/** Electron-only shell. Entering Expert mode requires no account/session action. */
export function DesktopProductShell({ flashClient, expertOwner }: DesktopProductShellProps) {
  const [mode, setMode] = useState<DesktopProductMode>('pt');
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
