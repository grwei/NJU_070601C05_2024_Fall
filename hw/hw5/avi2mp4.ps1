# 定义源目录路径
$sourceDirectory = "./fig/movie"

# 检查源目录是否存在
if (-Not (Test-Path $sourceDirectory)) {
    Write-Error "Source directory does not exist: $sourceDirectory"
    exit 1
}

# 获取所有匹配的 AVI 文件
$aviFiles = Get-ChildItem -Path $sourceDirectory -Filter "mov_*.avi" -File

# 遍历每个 AVI 文件并进行转换
foreach ($file in $aviFiles) {
    # 构建输出文件名（替换扩展名为 .mp4）
    $outputFilePath = [System.IO.Path]::ChangeExtension($file.FullName, ".mp4")

    # 使用 ffmpeg 进行无损转换
    & ffmpeg -i $file.FullName -c:v libx264 -crf 0 -preset slow -c:a aac -b:a 320k $outputFilePath

    # 检查转换是否成功
    if ($LASTEXITCODE -eq 0) {
        Write-Output "Successfully converted $($file.Name) to $([System.IO.Path]::GetFileName($outputFilePath))"
    } else {
        Write-Error "Failed to convert $($file.Name)"
    }
}
