$ErrorActionPreference = 'Stop'

param(
    [switch]$Hold,
    [switch]$NoHold
)

function Wait-IfNeeded {
    if ($NoHold) {
        return
    }

    if ($Hold) {
        Read-Host "按回车退出"
        return
    }

    if ($Host.Name -match 'ConsoleHost') {
        Read-Host "按回车退出"
    }
}

# One-click PETSc TS benchmark sweep for atmospheric_chemistry chamber S1.
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Resolve-Path (Join-Path $scriptDir "../..")
$actionsDir = Join-Path $repoRoot "modules/atmospheric_chemistry/test/tests/actions"

$candidates = @(
    (Join-Path $repoRoot "atmospheric_chemistry-opt.exe"),
    (Join-Path $repoRoot "atmospheric_chemistry-opt"),
    (Join-Path $repoRoot "modules/atmospheric_chemistry/atmospheric_chemistry-opt.exe"),
    (Join-Path $repoRoot "modules/atmospheric_chemistry/atmospheric_chemistry-opt")
)

$app = $null
foreach ($c in $candidates) {
    if (Test-Path $c) {
        $app = $c
        break
    }
}

if (-not $app) {
    Write-Error "atmospheric_chemistry-opt not found. Checked: $($candidates -join ', ')"
    Wait-IfNeeded
    exit 1
}

function Run-Case {
    param(
        [Parameter(Mandatory = $true)][string]$InputFile,
        [Parameter(Mandatory = $true)][string]$LogFile
    )

    Write-Host "===== Running $InputFile ====="
    & $app -i $InputFile -PerformanceLog 2>&1 | Tee-Object -FilePath $LogFile

    Write-Host "----- Key lines ($LogFile) -----"
    Select-String -Path $LogFile -Pattern "TOTAL RUN TIME IS|RUNNING F0AM|SAVED AS|TS" | ForEach-Object {
        $_.Line
    }
}

Set-Location $actionsDir
Run-Case -InputFile "vs_F0AM_chamber_S1_box_ts_bdf_rtol1e-2.i" -LogFile "ts_bdf_rtol1e-2.log"
Run-Case -InputFile "vs_F0AM_chamber_S1_box_ts_bdf_rtol1e-3.i" -LogFile "ts_bdf_rtol1e-3.log"
Run-Case -InputFile "vs_F0AM_chamber_S1_box_ts_arkimex_rtol1e-2.i" -LogFile "ts_arkimex_rtol1e-2.log"

Write-Host "===== Sweep complete ====="
Wait-IfNeeded
