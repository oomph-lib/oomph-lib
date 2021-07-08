<?php

echo '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
<head>
<center>
<img src="../../../doc/figures/oomph_logo.png" alt="oomph-lib logo" border=0></a>
</center>
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
$algebra=$_POST['algebra'];


$body="Name: $name\n\n";
$body.="Email: $email\n\n";
$body.="Affiliation: $affiliation\n\n";
$body.="Want to be put in mail list?: $maillist\n\n";
$body.="Why oomph?: $why\n\n";
$body.="First download: $firstdownload\n\n";
$body.="Algebra: $algebra\n\n";

$correct_algebra="seven";
if ( trim($algebra) != trim($correct_algebra) ) {
 die('Wrong spam filter answer...');
}
echo "Correct spam filter answer!";



$subject="oomph-lib 1.0.* download \n" ;

$headers="From: <$alias$email>\n"; 
$headers.="To: <oomphlib@maths.manchester.ac.uk>\n";


mail("", $subject, $body, $headers);
		
echo '
<body>






    Thank you! Now follow <B>
    <A HREF="../../../tar_file_directory">
    this link</A></B> to obtain the tar files. [Note for developers: This only works in the installed version.]

    <P>
    You may also want to download tar files that allow you install the third-party linear solvers
    trilinos, hypre and mumps as part of your installation. They are available here:
    <a href="../../../../oomph-lib_external_distfiles/">here</a>.

    <P>
 
    When you have downloaded the distribution (and possibly the third-party tar files) return to 
    the <A HREF="../../the_distribution/html/index.html#download">installation 
    instructions.</A>  


</body>
</html>





'

?>