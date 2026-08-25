$xml = [xml](Get-Content -Raw "D:\ONT\_xlsx_xl_worksheets_sheet1.xml")
$ns = New-Object System.Xml.XmlNamespaceManager($xml.NameTable)
$ns.AddNamespace("a","http://schemas.openxmlformats.org/spreadsheetml/2006/main")
$rows = $xml.SelectNodes("//a:sheetData/a:row", $ns)
$cols = @("Sample_ID","Group1","Group2","Gender","Age","Gravida","Loss","Para","Gestational_Week","Group3","CRL","Complications","Group4")
$data = New-Object System.Collections.Generic.List[object]
foreach ($row in $rows) {
  $vals = @{}
  foreach ($c in $row.SelectNodes("a:c", $ns)) {
    $ref = $c.GetAttribute("r")
    $t = $c.GetAttribute("t")
    if ($t -eq "inlineStr") { $v = $c.SelectSingleNode("a:is/a:t", $ns).InnerText }
    else { $v = $c.SelectSingleNode("a:v", $ns).InnerText }
    if ($ref -match '^([A-Z]+)(\d+)$') {
      $colLetter = $Matches[1]
      $colIdx = 0
      foreach ($ch in $colLetter.ToCharArray()) { $colIdx = $colIdx * 26 + ([int]$ch - 64) }
      $vals[$colIdx] = $v
    }
  }
  $line = New-Object System.Collections.Generic.List[string]
  for ($j=1; $j -le $cols.Count; $j++) {
    if ($vals.ContainsKey($j)) { $line.Add([string]$vals[$j]) } else { $line.Add("") }
  }
  $data.Add(($line -join "`t"))
}
[System.IO.File]::WriteAllLines("D:\ONT\clinical_649.tsv", $data, (New-Object System.Text.UTF8Encoding($false)))
"rows_written=" + $data.Count
# summary counts
$hdr = $data[0]
$body = $data[1..($data.Count-1)]
foreach ($ci in @(1,2,3,11,12)) {
  $name = $cols[$ci-1]
  $groups = $body | ForEach-Object { $_.Split("`t")[$ci-1] } | Group-Object | Sort-Object Count -Descending
  $s = ($groups | ForEach-Object { "$($_.Name)=$($_.Count)" }) -join ", "
  "[$name] $s"
}
