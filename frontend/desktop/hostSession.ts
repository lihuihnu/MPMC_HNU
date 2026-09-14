import { spawn, type ChildProcessWithoutNullStreams } from 'node:child_process';
import { randomBytes } from 'node:crypto';
import { isAbsolute } from 'node:path';
import { createInterface } from 'node:readline';

const HOST_START_TIMEOUT_MS = 20_000;
const HOST_STOP_TIMEOUT_MS = 12_000;
const MAX_HOST_LOG_LINE = 64 * 1024;

export interface PtHostConnection {
  readonly baseUrl: string;
  readonly bearerToken: string;
  readonly configuredBackends: number;
}

interface PtReadyEvent {
  event: 'pt_process_ready';
  selected_port: number;
  configured_backends: number;
}

export function parsePtReadyEvent(line: string): PtReadyEvent | null {
  if (line.length === 0 || line.length > MAX_HOST_LOG_LINE) {
    return null;
  }
  let value: unknown;
  try {
    value = JSON.parse(line);
  } catch {
    return null;
  }
  if (typeof value !== 'object' || value === null) {
    return null;
  }
  const event = value as Partial<PtReadyEvent>;
  if (
    event.event !== 'pt_process_ready' ||
    !Number.isInteger(event.selected_port) ||
    event.selected_port === undefined ||
    event.selected_port < 1 ||
    event.selected_port > 65_535 ||
    !Number.isInteger(event.configured_backends) ||
    event.configured_backends !== 3
  ) {
    return null;
  }
  return {
    event: 'pt_process_ready',
    selected_port: event.selected_port,
    configured_backends: event.configured_backends,
  } as PtReadyEvent;
}

export class PtHostSession {
  private child: ChildProcessWithoutNullStreams | null = null;
  private connection: PtHostConnection | null = null;
  private starting: Promise<PtHostConnection> | null = null;

  constructor(
    private readonly binaryPath: string,
    private readonly log: (line: string) => void = () => {},
  ) {
    if (!isAbsolute(binaryPath)) {
      throw new Error('The PT desktop host path must be absolute.');
    }
  }

  start(): Promise<PtHostConnection> {
    if (this.connection !== null) {
      return Promise.resolve(this.connection);
    }
    if (this.starting !== null) {
      return this.starting;
    }
    this.starting = this.startChild().finally(() => {
      this.starting = null;
    });
    return this.starting;
  }

  private startChild(): Promise<PtHostConnection> {
    const token = randomBytes(32).toString('base64url');
    const child = spawn(this.binaryPath, ['--desktop-session-token-stdin'], {
      detached: false,
      shell: false,
      stdio: ['pipe', 'pipe', 'pipe'],
      windowsHide: true,
    });
    this.child = child;
    child.stdin.setDefaultEncoding('utf8');

    const stdout = createInterface({ input: child.stdout, crlfDelay: Infinity });
    const stderr = createInterface({ input: child.stderr, crlfDelay: Infinity });
    stdout.on('line', (line) => {
      if (line.length <= MAX_HOST_LOG_LINE) {
        this.log(line);
      }
    });
    stderr.on('line', (line) => {
      if (line.length <= MAX_HOST_LOG_LINE) {
        this.log(line);
      }
    });
    child.once('exit', () => {
      if (this.child === child) {
        this.child = null;
        this.connection = null;
      }
    });

    return new Promise<PtHostConnection>((resolve, reject) => {
      let settled = false;
      const finishError = (error: Error) => {
        if (settled) {
          return;
        }
        settled = true;
        clearTimeout(timeout);
        child.stdin.destroy();
        child.kill();
        reject(error);
      };
      const timeout = setTimeout(() => {
        finishError(new Error('The PT desktop host did not become ready in time.'));
      }, HOST_START_TIMEOUT_MS);

      child.once('error', (error) => {
        finishError(new Error(`The PT desktop host could not start: ${error.message}`));
      });
      child.stdin.once('error', (error) => {
        finishError(
          new Error(`The PT desktop session token could not be delivered: ${error.message}`),
        );
      });
      child.once('exit', (code, signal) => {
        finishError(
          new Error(
            `The PT desktop host exited before readiness (code ${String(code)}, signal ${String(signal)}).`,
          ),
        );
      });
      stdout.on('line', (line) => {
        const ready = parsePtReadyEvent(line);
        if (settled || ready === null) {
          return;
        }
        settled = true;
        clearTimeout(timeout);
        const connection: PtHostConnection = {
          baseUrl: `http://127.0.0.1:${ready.selected_port}`,
          bearerToken: token,
          configuredBackends: ready.configured_backends,
        };
        this.connection = connection;
        resolve(connection);
      });
      child.stdin.write(`${token}\n`);
    });
  }

  async stop(): Promise<void> {
    const child = this.child;
    this.connection = null;
    if (child === null) {
      return;
    }
    this.child = null;
    if (child.exitCode !== null || child.signalCode !== null) {
      return;
    }
    const exited = new Promise<void>((resolve) => {
      child.once('exit', () => resolve());
    });
    if (!child.killed && child.stdin.writable) {
      child.stdin.end('shutdown\n');
    }
    let stopTimer: ReturnType<typeof setTimeout> | undefined;
    const graceful = await Promise.race([
      exited.then(() => true),
      new Promise<false>((resolve) => {
        stopTimer = setTimeout(() => resolve(false), HOST_STOP_TIMEOUT_MS);
      }),
    ]);
    if (stopTimer !== undefined) {
      clearTimeout(stopTimer);
    }
    if (!graceful && !child.killed) {
      child.kill();
      await exited;
    }
  }
}
