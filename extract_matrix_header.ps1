$path = "E:\甲基化数据矩阵\EWAS_INPUT_NO_HEADER.txt"
$sr = New-Object System.IO.StreamReader($path, [System.Text.Encoding]::UTF8)
$rows = @()
for ($i=0; $i -lt 71; $i++) {
  if ($sr.EndOfStream) { break }
  $rows += ,($sr.ReadLine().Split("`t"))
}
$sr.Close()
"rows_read = " + $rows.Count
$nSamples = $rows[0].Count - 1
"nSamples = " + $nSamples
$labels = @()
foreach ($r in $rows) { $labels += $r[0].Substring(1) }
$out = New-Object System.Collections.Generic.List[string]
$out.Add(($labels -join "`t"))
for ($j=0; $j -lt $nSamples; $j++) {
  $vals = @()
  foreach ($r in $rows) { $vals += $r[$j+1] }
  $out.Add(($vals -join "`t"))
}
[System.IO.File]::WriteAllLines("D:\ONT\matrix_covariates.tsv", $out, (New-Object System.Text.UTF8Encoding($false)))
$fidLine = $rows[1]
$fidList = @($fidLine[1..($fidLine.Count-1)])
[System.IO.File]::WriteAllLines("D:\ONT\matrix_fid.txt", $fidList, (New-Object System.Text.UTF8Encoding($false)))
"wrote covariates rows = " + $out.Count
"wrote fid list count = " + $fidList.Count
