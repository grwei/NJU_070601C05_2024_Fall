$overwriteAll = $false

# 获取所有匹配 ./hwN/hwN.html 模式的文件
$files = Get-ChildItem -Path .\hw -Recurse -Filter "hw*.html" | Where-Object {
    $_.FullName -match 'hw(\d+)\\hw\1\.html$'
}

foreach ($file in $files) {
    $sourceFile = $file.FullName
    $destinationFile = Join-Path -Path "./docs" -ChildPath (Split-Path -Leaf $sourceFile)

    if (Test-Path $destinationFile) {
        if (-not $overwriteAll) {
            Write-Host "The file $destinationFile already exists. What would you like to do?"
            Write-Host "[Y] Yes, overwrite the current file."
            Write-Host "[A] Yes, overwrite all files without asking."
            Write-Host "[N] No, skip this file."

            $response = Read-Host "Please select an option [Y/A/N]"

            switch ($response.ToUpper()) {
                'Y' {
                    Copy-Item -Path $sourceFile -Destination $destinationFile -Force
                    Write-Host "$destinationFile has been overwritten."
                }
                'A' {
                    $overwriteAll = $true
                    Copy-Item -Path $sourceFile -Destination $destinationFile -Force
                    Write-Host "$destinationFile has been overwritten."
                }
                'N' {
                    Write-Host "$destinationFile will be skipped."
                    continue
                }
                default {
                    Write-Host "Invalid response. Skipping $destinationFile."
                    continue
                }
            }
        }
        else {
            Copy-Item -Path $sourceFile -Destination $destinationFile -Force
            Write-Host "$destinationFile has been overwritten."
        }
    }
    else {
        Copy-Item -Path $sourceFile -Destination $destinationFile
        Write-Host "$destinationFile has been copied."
    }
}
