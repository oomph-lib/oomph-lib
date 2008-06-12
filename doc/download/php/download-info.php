<?php

echo '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
<link href="../../oomph_header.css" rel="stylesheet" type="text/css"> 
<a href="http://www.oomph-lib.org"><img src="../../../doc/figures/oomph.png" alt="oomph-lib logo" border=0></a>
<HR>
<title>Oomph-lib download page</title>
</head>
';
$name=$_POST['name'];
$affiliation=$_POST['affiliation'];
$email=$_POST['email'];
$maillist=$_POST['maillist'];
$why=$_POST['why'];
$firstdownload=$_POST['firstdownload'];


$body="Name: $name\n\n";
$body.="Email: $email\n\n";
$body.="Affiliation: $affiliation\n\n";
$body.="Want to be put in mail list?: $maillist\n\n";
$body.="Why oomph?: $why\n\n";
$body.="First download: $firstdownload\n\n";



$subject="oomph download \n" ;

$headers="From: <$alias$email>\n"; 
$headers.="To: <mheil@maths.man.ac.uk>\n";


mail("", $subject, $body, $headers);
		
echo '
<body>


<p>Choose the distribution:</p>


<CENTER>
<TABLE BORDER=1>
<TR>
<TD>
<B>Tar file</B>
</TD>
<TD>
<B>Content</B>
</TD>
<TD>
<B>Size of the gzipped tar file</B>
</TD>
<TD>
<B>Size of the built distribution</B>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../tar_files/html/index.html">
oomph-lib-0.0.tar.gz
</A>
</TD>
<TD>
Everything:
- sources
- demo driver codes
- documentation
</TD>
<TD>
10 gazillion MB
</TD>
<TD>
20 gazillion MB
</TD>
</TR>
<TR>
<TD>
<A HREF="../../tar_files/html/index.html">
oomph-lib-no-doc-0.0.tar.gz
</A>
</TD>
<TD>
Everything apart from the documentation:
- sources
- demo driver codes
</TD>
<TD>
1 gazillion MB
</TD>
<TD>
2 gazillion MB
</TD>
</TR>
<TR>
<TD>
<A HREF="../../tar_files/html/index.html">
oomph-lib-no-doc-no-demo-0.0.tar.gz
</A>
</TD>
<TD>
Only the library sources:
- sources
</TD>
<TD>
0.1 gazillion MB
</TD>
<TD>
0.2 gazillion MB
</TD>
</TR>
</TABLE>
</CENTER>

<P>

When you have downloaded the distribution return to 
the <A HREF="../../the_distribution/html/index.html#download">installation 
instructions.</A>

</P>

</body>
</html>




'

?>