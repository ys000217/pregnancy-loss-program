$xml = [xml](Get-Content -Raw "D:\ONT\_xlsx_xl_worksheets_sheet1.xml")
$ns = New-Object System.Xml.XmlNamespaceManager($xml.NameTable)
$ns.AddNamespace("a","http://schemas.openxmlformats.org/spreadsheetml/2006/main")
$rows = $xml.SelectNodes("//a:sheetData/a:row", $ns)
$n = [Math]::Min(8, $rows.Count)
for ($i=0; $i -lt $n; $i++) {
  $row = $rows[$i]
  $cells = @()
  foreach ($c in $row.SelectNodes("a:c", $ns)) {
    $ref = $c.GetAttribute("r")
    $t = $c.GetAttribute("t")
    if ($t -eq "inlineStr") { $v = $c.SelectSingleNode("a:is/a:t", $ns).InnerText }
    else { $v = $c.SelectSingleNode("a:v", $ns).InnerText }
    $cells += "${ref}=${v}"
  }
  "ROW$($i+1): " + ($cells -join " | ")
}
