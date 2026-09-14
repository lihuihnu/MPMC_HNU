[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('payload', 'msi')]
    [string] $Mode,

    [Parameter(Mandatory = $false)]
    [string] $PackageDir,

    [Parameter(Mandatory = $false)]
    [string] $MsiPath,

    [Parameter(Mandatory = $true)]
    [string] $EvidencePath
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$codeSigningOid = '1.3.6.1.5.5.7.3.3'
$requiredEnvironment = @(
    'MPMC_WINDOWS_CODESIGN_PFX_PATH',
    'MPMC_WINDOWS_CODESIGN_PFX_PASSWORD',
    'MPMC_WINDOWS_CODESIGN_TIMESTAMP_URL',
    'MPMC_WINDOWS_CODESIGN_EXPECTED_SUBJECT'
)

foreach ($name in $requiredEnvironment) {
    $value = [Environment]::GetEnvironmentVariable($name)
    if ([string]::IsNullOrWhiteSpace($value)) {
        throw "Required Windows release signing setting is missing: $name"
    }
}

$pfxPath = [Environment]::GetEnvironmentVariable('MPMC_WINDOWS_CODESIGN_PFX_PATH')
$pfxPassword = [Environment]::GetEnvironmentVariable('MPMC_WINDOWS_CODESIGN_PFX_PASSWORD')
$timestampUrl = [Environment]::GetEnvironmentVariable('MPMC_WINDOWS_CODESIGN_TIMESTAMP_URL')
$expectedSubject = [Environment]::GetEnvironmentVariable('MPMC_WINDOWS_CODESIGN_EXPECTED_SUBJECT')

if (-not (Test-Path $pfxPath -PathType Leaf)) {
    throw 'The configured Authenticode PFX file does not exist'
}
try {
    $timestampUri = [Uri] $timestampUrl
} catch {
    throw 'The Authenticode RFC 3161 timestamp URL is invalid'
}
if ($timestampUri.Scheme -notin @('http', 'https')) {
    throw 'The Authenticode timestamp URL must use HTTP or HTTPS'
}

$signTool = Get-ChildItem -Path "${env:ProgramFiles(x86)}\Windows Kits\10\bin\*\x64\signtool.exe" |
    Sort-Object FullName -Descending |
    Select-Object -First 1
if ($null -eq $signTool) {
    throw 'Windows SDK SignTool is unavailable'
}

$securePassword = ConvertTo-SecureString $pfxPassword -AsPlainText -Force
$imported = @(
    Import-PfxCertificate -FilePath $pfxPath -CertStoreLocation 'Cert:\CurrentUser\My' -Password $securePassword -Exportable:$false
)

try {
    $certificate = $imported |
        Where-Object {
            $_.HasPrivateKey -and
            @($_.EnhancedKeyUsageList | ForEach-Object { $_.ObjectId.Value }) -contains $codeSigningOid
        } |
        Select-Object -First 1
    if ($null -eq $certificate) {
        throw 'The supplied PFX contains no private-key code-signing certificate'
    }
    if ($certificate.Subject -ne $expectedSubject) {
        throw "Code-signing certificate subject mismatch: $($certificate.Subject)"
    }
    $now = Get-Date
    if ($certificate.NotBefore -gt $now -or $certificate.NotAfter -le $now) {
        throw 'The supplied code-signing certificate is not currently valid'
    }

    $targets = @()
    if ($Mode -eq 'payload') {
        if ([string]::IsNullOrWhiteSpace($PackageDir)) {
            throw 'PackageDir is required for payload signing'
        }
        $packageRoot = (Resolve-Path $PackageDir).Path
        $targets = @(
            @{
                Role = 'desktop-executable'
                Path = Join-Path $packageRoot 'MPMC-PT-Desktop.exe'
                Relative = 'MPMC-PT-Desktop.exe'
            },
            @{
                Role = 'native-host'
                Path = Join-Path $packageRoot 'resources\desktop-native\bin\mpmc_pt_service_host.exe'
                Relative = 'resources/desktop-native/bin/mpmc_pt_service_host.exe'
            }
        )
    } else {
        if ([string]::IsNullOrWhiteSpace($MsiPath)) {
            throw 'MsiPath is required for MSI signing'
        }
        $resolvedMsi = (Resolve-Path $MsiPath).Path
        $targets = @(
            @{
                Role = 'msi'
                Path = $resolvedMsi
                Relative = [IO.Path]::GetFileName($resolvedMsi)
            }
        )
    }

    foreach ($target in $targets) {
        if (-not (Test-Path $target.Path -PathType Leaf)) {
            throw "Required Authenticode target is missing: $($target.Path)"
        }
    }

    $newEntries = @()
    foreach ($target in $targets) {
        $signArguments = @(
            'sign',
            '/sha1', $certificate.Thumbprint,
            '/fd', 'SHA256',
            '/tr', $timestampUrl,
            '/td', 'SHA256',
            '/v',
            $target.Path
        )
        & $signTool.FullName @signArguments
        if ($LASTEXITCODE -ne 0) {
            throw "SignTool failed to sign $($target.Role) with exit code $LASTEXITCODE"
        }

        & $signTool.FullName verify /pa /all /v $target.Path
        if ($LASTEXITCODE -ne 0) {
            throw "SignTool verification failed for $($target.Role) with exit code $LASTEXITCODE"
        }

        $signature = Get-AuthenticodeSignature -FilePath $target.Path
        if ($signature.Status -ne 'Valid') {
            throw "Authenticode status is not Valid for $($target.Role): $($signature.Status)"
        }
        if ($null -eq $signature.SignerCertificate -or
            $signature.SignerCertificate.Subject -ne $expectedSubject) {
            throw "Authenticode signer subject mismatch for $($target.Role)"
        }
        if ($null -eq $signature.TimeStamperCertificate) {
            throw "RFC 3161 timestamp evidence is missing for $($target.Role)"
        }

        $newEntries += [ordered]@{
            role = $target.Role
            file = $target.Relative
            sha256 = (Get-FileHash -Algorithm SHA256 -Path $target.Path).Hash.ToLowerInvariant()
            status = 'Valid'
            signer_subject = $signature.SignerCertificate.Subject
            signer_thumbprint = $signature.SignerCertificate.Thumbprint.ToLowerInvariant()
            timestamped = $true
            timestamp_subject = $signature.TimeStamperCertificate.Subject
        }
    }

    $evidenceFile = [IO.Path]::GetFullPath($EvidencePath)
    $existingEntries = @()
    if (Test-Path $evidenceFile -PathType Leaf) {
        $existing = Get-Content $evidenceFile -Raw | ConvertFrom-Json
        if ($existing.convention -ne 'MPMC/PT/windows-authenticode-evidence/v1') {
            throw 'Existing Authenticode evidence convention changed'
        }
        if ($existing.signer_subject -ne $expectedSubject) {
            throw 'Existing Authenticode evidence signer changed'
        }
        $existingEntries = @($existing.files)
    }

    $rolesBeingReplaced = @($newEntries | ForEach-Object { $_.role })
    $mergedEntries = @(
        $existingEntries | Where-Object { $_.role -notin $rolesBeingReplaced }
    ) + $newEntries

    $evidence = [ordered]@{
        convention = 'MPMC/PT/windows-authenticode-evidence/v1'
        signer_subject = $expectedSubject
        signer_thumbprint = $certificate.Thumbprint.ToLowerInvariant()
        file_digest = 'SHA256'
        timestamp_protocol = 'RFC3161'
        timestamp_digest = 'SHA256'
        files = @($mergedEntries | Sort-Object role)
    }
    $parent = Split-Path -Parent $evidenceFile
    if (-not [string]::IsNullOrWhiteSpace($parent)) {
        New-Item -ItemType Directory -Path $parent -Force | Out-Null
    }
    $evidence | ConvertTo-Json -Depth 8 | Set-Content -Path $evidenceFile -Encoding utf8
    Write-Host "WINDOWS_AUTHENTICODE_$($Mode.ToUpperInvariant())_OK"
} finally {
    foreach ($item in $imported) {
        if ($null -ne $item.Thumbprint) {
            Remove-Item "Cert:\CurrentUser\My\$($item.Thumbprint)" -Force -ErrorAction SilentlyContinue
        }
    }
}
