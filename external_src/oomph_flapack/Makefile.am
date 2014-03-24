###################################################
#   FORTRAN LAPACK sources
###################################################

# Define the sources
sources = dummy.cc \
cbdsqr.f  cpptrf.f    dlaqps.f  sbdsdc.f    sorm2l.f  zhetrf.f \
cgbbrd.f  cpptri.f    dlaqr0.f  sbdsqr.f    sorm2r.f  zhetri.f \
cgbcon.f  cpptrs.f    dlaqr1.f  scsum1.f    sormbr.f  zhetrs.f \
cgbequ.f  cptcon.f    dlaqr2.f  sdisna.f    sormhr.f  zhgeqz.f \
cgbrfs.f  cpteqr.f    dlaqr3.f  sgbbrd.f    sorml2.f  zhpcon.f \
cgbsv.f   cptrfs.f    dlaqr4.f  sgbcon.f    sormlq.f  zhpevd.f \
cgbsvx.f  cptsv.f     dlaqr5.f  sgbequ.f    sormql.f  zhpev.f \
cgbtf2.f  cptsvx.f    dlaqsb.f  sgbrfs.f    sormqr.f  zhpevx.f \
cgbtrf.f  cpttrf.f    dlaqsp.f  sgbsv.f     sormr2.f  zhpgst.f \
cgbtrs.f  cpttrs.f    dlaqsy.f  sgbsvx.f    sormr3.f  zhpgvd.f \
cgebak.f  cptts2.f    dlaqtr.f  sgbtf2.f    sormrq.f  zhpgv.f \
cgebal.f  crot.f      dlar1v.f  sgbtrf.f    sormrz.f  zhpgvx.f \
cgebd2.f  cspcon.f    dlar2v.f  sgbtrs.f    sormtr.f  zhprfs.f \
cgebrd.f  cspmv.f     dlaran.f  sgebak.f    spbcon.f  zhpsv.f \
cgecon.f  cspr.f      dlarfb.f  sgebal.f    spbequ.f  zhpsvx.f \
cgeequ.f  csprfs.f    dlarf.f   sgebd2.f    spbrfs.f  zhptrd.f \
cgees.f   cspsv.f     dlarfg.f  sgebrd.f    spbstf.f  zhptrf.f \
cgeesx.f  cspsvx.f    dlarft.f  sgecon.f    spbsv.f   zhptri.f \
cgeev.f   csptrf.f    dlarfx.f  sgeequ.f    spbsvx.f  zhptrs.f \
cgeevx.f  csptri.f    dlargv.f  sgees.f     spbtf2.f  zhsein.f \
cgegs.f   csptrs.f    dlarnd.f  sgeesx.f    spbtrf.f  zhseqr.f \
cgegv.f   csrot.f     dlarnv.f  sgeev.f     spbtrs.f  zlabrd.f \
cgehd2.f  csrscl.f    dlarra.f  sgeevx.f    spocon.f  zlacgv.f \
cgehrd.f  cstedc.f    dlarrb.f  sgegs.f     spoequ.f  zlacn2.f \
cgelq2.f  cstegr.f    dlarrc.f  sgegv.f     sporfs.f  zlacon.f \
cgelqf.f  cstein.f    dlarrd.f  sgehd2.f    sposv.f   zlacp2.f \
cgelsd.f  cstemr.f    dlarre.f  sgehrd.f    sposvx.f  zlacpy.f \
cgels.f   csteqr.f    dlarrf.f  sgelq2.f    spotf2.f  zlacrm.f \
cgelss.f  csycon.f    dlarrj.f  sgelqf.f    spotrf.f  zlacrt.f \
cgelsx.f  csymv.f     dlarrk.f  sgelsd.f    spotri.f  zladiv.f \
cgelsy.f  csyr.f      dlarrr.f  sgels.f     spotrs.f  zlaed0.f \
cgeql2.f  csyrfs.f    dlarrv.f  sgelss.f    sppcon.f  zlaed7.f \
cgeqlf.f  csysv.f     dlartg.f  sgelsx.f    sppequ.f  zlaed8.f \
cgeqp3.f  csysvx.f    dlartv.f  sgelsy.f    spprfs.f  zlaein.f \
cgeqpf.f  csytf2.f    dlaruv.f  sgeql2.f    sppsv.f   zlaesy.f \
cgeqr2.f  csytrf.f    dlarzb.f  sgeqlf.f    sppsvx.f  zlaev2.f \
cgeqrf.f  csytri.f    dlarz.f   sgeqp3.f    spptrf.f  zlag2c.f \
cgerfs.f  csytrs.f    dlarzt.f  sgeqpf.f    spptri.f  zlags2.f \
cgerq2.f  ctbcon.f    dlas2.f   sgeqr2.f    spptrs.f  zlagtm.f \
cgerqf.f  ctbrfs.f    dlascl.f  sgeqrf.f    sptcon.f  zlahef.f \
cgesc2.f  ctbtrs.f    dlasd0.f  sgerfs.f    spteqr.f  zlahqr.f \
cgesdd.f  ctgevc.f    dlasd1.f  sgerq2.f    sptrfs.f  zlahr2.f \
cgesvd.f  ctgex2.f    dlasd2.f  sgerqf.f    sptsv.f   zlahrd.f \
cgesv.f   ctgexc.f    dlasd3.f  sgesc2.f    sptsvx.f  zlaic1.f \
cgesvx.f  ctgsen.f    dlasd4.f  sgesdd.f    spttrf.f  zlals0.f \
cgetc2.f  ctgsja.f    dlasd5.f  sgesvd.f    spttrs.f  zlalsa.f \
cgetf2.f  ctgsna.f    dlasd6.f  sgesv.f     sptts2.f  zlalsd.f \
cgetrf.f  ctgsy2.f    dlasd7.f  sgesvx.f    srscl.f   zlangb.f \
cgetri.f  ctgsyl.f    dlasd8.f  sgetc2.f    ssbevd.f  zlange.f \
cgetrs.f  ctpcon.f    dlasda.f  sgetf2.f    ssbev.f   zlangt.f \
cggbak.f  ctprfs.f    dlasdq.f  sgetrf.f    ssbevx.f  zlanhb.f \
cggbal.f  ctptri.f    dlasdt.f  sgetri.f    ssbgst.f  zlanhe.f \
cgges.f   ctptrs.f    dlaset.f  sgetrs.f    ssbgvd.f  zlanhp.f \
cggesx.f  ctrcon.f    dlasq1.f  sggbak.f    ssbgv.f   zlanhs.f \
cggev.f   ctrevc.f    dlasq2.f  sggbal.f    ssbgvx.f  zlanht.f \
cggevx.f  ctrexc.f    dlasq3.f  sgges.f     ssbtrd.f  zlansb.f \
cggglm.f  ctrrfs.f    dlasq4.f  sggesx.f    sspcon.f  zlansp.f \
cgghrd.f  ctrsen.f    dlasq5.f  sggev.f     sspevd.f  zlansy.f \
cgglse.f  ctrsna.f    dlasq6.f  sggevx.f    sspev.f   zlantb.f \
cggqrf.f  ctrsyl.f    dlasr.f   sggglm.f    sspevx.f  zlantp.f \
cggrqf.f  ctrti2.f    dlasrt.f  sgghrd.f    sspgst.f  zlantr.f \
cggsvd.f  ctrtri.f    dlassq.f  sgglse.f    sspgvd.f  zlapll.f \
cggsvp.f  ctrtrs.f    dlasv2.f  sggqrf.f    sspgv.f   zlapmt.f \
cgtcon.f  ctzrqf.f    dlaswp.f  sggrqf.f    sspgvx.f  zlaqgb.f \
cgtrfs.f  ctzrzf.f    dlasy2.f  sggsvd.f    ssprfs.f  zlaqge.f \
cgtsv.f   cung2l.f    dlasyf.f  sggsvp.f    sspsv.f   zlaqhb.f \
cgtsvx.f  cung2r.f    dlatbs.f  sgtcon.f    sspsvx.f  zlaqhe.f \
cgttrf.f  cungbr.f    dlatdf.f  sgtrfs.f    ssptrd.f  zlaqhp.f \
cgttrs.f  cunghr.f    dlatps.f  sgtsv.f     ssptrf.f  zlaqp2.f \
cgtts2.f  cungl2.f    dlatrd.f  sgtsvx.f    ssptri.f  zlaqps.f \
chbevd.f  cunglq.f    dlatrs.f  sgttrf.f    ssptrs.f  zlaqr0.f \
chbev.f   cungql.f    dlatrz.f  sgttrs.f    sstebz.f  zlaqr1.f \
chbevx.f  cungqr.f    dlatzm.f  sgtts2.f    sstedc.f  zlaqr2.f \
chbgst.f  cungr2.f    dlauu2.f  shgeqz.f    sstegr.f  zlaqr3.f \
chbgvd.f  cungrq.f    dlauum.f  shsein.f    sstein.f  zlaqr4.f \
chbgv.f   cungtr.f    dlazq3.f  shseqr.f    sstemr.f  zlaqr5.f \
chbgvx.f  cunm2l.f    dlazq4.f  sisnan.f    ssteqr.f  zlaqsb.f \
chbtrd.f  cunm2r.f    dopgtr.f  slabad.f    ssterf.f  zlaqsp.f \
checon.f  cunmbr.f    dopmtr.f  slabrd.f    sstevd.f  zlaqsy.f \
cheevd.f  cunmhr.f    dorg2l.f  slacn2.f    sstev.f   zlar1v.f \
cheev.f   cunml2.f    dorg2r.f  slacon.f    sstevr.f  zlar2v.f \
cheevr.f  cunmlq.f    dorgbr.f  slacpy.f    sstevx.f  zlarcm.f \
cheevx.f  cunmql.f    dorghr.f  sladiv.f    ssycon.f  zlarfb.f \
chegs2.f  cunmqr.f    dorgl2.f  slae2.f     ssyevd.f  zlarf.f \
chegst.f  cunmr2.f    dorglq.f  slaebz.f    ssyev.f   zlarfg.f \
chegvd.f  cunmr3.f    dorgql.f  slaed0.f    ssyevr.f  zlarft.f \
chegv.f   cunmrq.f    dorgqr.f  slaed1.f    ssyevx.f  zlarfx.f \
chegvx.f  cunmrz.f    dorgr2.f  slaed2.f    ssygs2.f  zlargv.f \
cherfs.f  cunmtr.f    dorgrq.f  slaed3.f    ssygst.f  zlarnv.f \
chesv.f   cupgtr.f    dorgtr.f  slaed4.f    ssygvd.f  zlarrv.f \
chesvx.f  cupmtr.f    dorm2l.f  slaed5.f    ssygv.f   zlartg.f \
chetd2.f  dbdsdc.f    dorm2r.f  slaed6.f    ssygvx.f  zlartv.f \
chetf2.f  dbdsqr.f    dormbr.f  slaed7.f    ssyrfs.f  zlarzb.f \
chetrd.f  ddisna.f    dormhr.f  slaed8.f    ssysv.f   zlarz.f \
chetrf.f  dgbbrd.f    dorml2.f  slaed9.f    ssysvx.f  zlarzt.f \
chetri.f  dgbcon.f    dormlq.f  slaeda.f    ssytd2.f  zlascl.f \
chetrs.f  dgbequ.f    dormql.f  slaein.f    ssytf2.f  zlaset.f \
chgeqz.f  dgbrfs.f    dormqr.f  slaev2.f    ssytrd.f  zlasr.f \
chpcon.f  dgbsv.f     dormr2.f  slaexc.f    ssytrf.f  zlassq.f \
chpevd.f  dgbsvx.f    dormr3.f  slag2d.f    ssytri.f  zlaswp.f \
chpev.f   dgbtf2.f    dormrq.f  slag2.f     ssytrs.f  zlasyf.f \
chpevx.f  dgbtrf.f    dormrz.f  slags2.f    stbcon.f  zlatbs.f \
chpgst.f  dgbtrs.f    dormtr.f  slagtf.f    stbrfs.f  zlatdf.f \
chpgvd.f  dgebak.f    dpbcon.f  slagtm.f    stbtrs.f  zlatps.f \
chpgv.f   dgebal.f    dpbequ.f  slagts.f    stgevc.f  zlatrd.f \
chpgvx.f  dgebd2.f    dpbrfs.f  slagv2.f    stgex2.f  zlatrs.f \
chprfs.f  dgebrd.f    dpbstf.f  slahqr.f    stgexc.f  zlatrz.f \
chpsv.f   dgecon.f    dpbsv.f   slahr2.f    stgsen.f  zlatzm.f \
chpsvx.f  dgeequ.f    dpbsvx.f  slahrd.f    stgsja.f  zlauu2.f \
chptrd.f  dgees.f     dpbtf2.f  slaic1.f    stgsna.f  zlauum.f \
chptrf.f  dgeesx.f    dpbtrf.f  slaisnan.f  stgsy2.f  zmach.f \
chptri.f  dgeev.f     dpbtrs.f  slaln2.f    stgsyl.f  zpbcon.f \
chptrs.f  dgeevx.f    dpocon.f  slals0.f    stpcon.f  zpbequ.f \
chsein.f  dgegs.f     dpoequ.f  slalsa.f    stprfs.f  zpbrfs.f \
chseqr.f  dgegv.f     dporfs.f  slalsd.f    stptri.f  zpbstf.f \
clabrd.f  dgehd2.f    dposv.f   slamch.f    stptrs.f  zpbsv.f \
clacgv.f  dgehrd.f    dposvx.f  slamrg.f    strcon.f  zpbsvx.f \
clacn2.f  dgelq2.f    dpotf2.f  slaneg.f    strevc.f  zpbtf2.f \
clacon.f  dgelqf.f    dpotrf.f  slangb.f    strexc.f  zpbtrf.f \
clacp2.f  dgelsd.f    dpotri.f  slange.f    strrfs.f  zpbtrs.f \
clacpy.f  dgels.f     dpotrs.f  slangt.f    strsen.f  zpocon.f \
clacrm.f  dgelss.f    dppcon.f  slanhs.f    strsna.f  zpoequ.f \
clacrt.f  dgelsx.f    dppequ.f  slansb.f    strsyl.f  zporfs.f \
cladiv.f  dgelsy.f    dpprfs.f  slansp.f    strti2.f  zposv.f \
claed0.f  dgeql2.f    dppsv.f   slanst.f    strtri.f  zposvx.f \
claed7.f  dgeqlf.f    dppsvx.f  slansy.f    strtrs.f  zpotf2.f \
claed8.f  dgeqp3.f    dpptrf.f  slantb.f    stzrqf.f  zpotrf.f \
claein.f  dgeqpf.f    dpptri.f  slantp.f    stzrzf.f  zpotri.f \
claesy.f  dgeqr2.f    dpptrs.f  slantr.f    xerbla.f  zpotrs.f \
claev2.f  dgeqrf.f    dptcon.f  slanv2.f    xlaenv.f  zppcon.f \
clag2z.f  dgerfs.f    dpteqr.f  slapll.f    zbdsqr.f  zppequ.f \
clags2.f  dgerq2.f    dptrfs.f  slapmt.f    zcgesv.f  zpprfs.f \
clagtm.f  dgerqf.f    dptsv.f   slaptm.f    zdrot.f   zppsv.f \
clahef.f  dgesc2.f    dptsvx.f  slapy2.f    zdrscl.f  zppsvx.f \
clahqr.f  dgesdd.f    dpttrf.f  slapy3.f    zdscal.f  zpptrf.f \
clahr2.f  dgesvd.f    dpttrs.f  slaqgb.f    zgbbrd.f  zpptri.f \
clahrd.f  dgesv.f     dptts2.f  slaqge.f    zgbcon.f  zpptrs.f \
claic1.f  dgesvx.f    drscl.f   slaqp2.f    zgbequ.f  zptcon.f \
clals0.f  dgetc2.f    dsbevd.f  slaqps.f    zgbrfs.f  zpteqr.f \
clalsa.f  dgetf2.f    dsbev.f   slaqr0.f    zgbsv.f   zptrfs.f \
clalsd.f  dgetrf.f    dsbevx.f  slaqr1.f    zgbsvx.f  zptsv.f \
clangb.f  dgetri.f    dsbgst.f  slaqr2.f    zgbtf2.f  zptsvx.f \
clange.f  dgetrs.f    dsbgvd.f  slaqr3.f    zgbtrf.f  zpttrf.f \
clangt.f  dggbak.f    dsbgv.f   slaqr4.f    zgbtrs.f  zpttrs.f \
clanhb.f  dggbal.f    dsbgvx.f  slaqr5.f    zgebak.f  zptts2.f \
clanhe.f  dgges.f     dsbtrd.f  slaqsb.f    zgebal.f  zrot.f \
clanhp.f  dggesx.f    dsgesv.f  slaqsp.f    zgebd2.f  zspcon.f \
clanhs.f  dggev.f     dspcon.f  slaqsy.f    zgebrd.f  zspmv.f \
clanht.f  dggevx.f    dspevd.f  slaqtr.f    zgecon.f  zspr.f \
clansb.f  dggglm.f    dspev.f   slar1v.f    zgeequ.f  zsprfs.f \
clansp.f  dgghrd.f    dspevx.f  slar2v.f    zgees.f   zspsv.f \
clansy.f  dgglse.f    dspgst.f  slaran.f    zgeesx.f  zspsvx.f \
clantb.f  dggqrf.f    dspgvd.f  slarfb.f    zgeev.f   zsptrf.f \
clantp.f  dggrqf.f    dspgv.f   slarf.f     zgeevx.f  zsptri.f \
clantr.f  dggsvd.f    dspgvx.f  slarfg.f    zgegs.f   zsptrs.f \
clapll.f  dggsvp.f    dsprfs.f  slarft.f    zgegv.f   zstedc.f \
clapmt.f  dgtcon.f    dspsv.f   slarfx.f    zgehd2.f  zstegr.f \
claqgb.f  dgtrfs.f    dspsvx.f  slargv.f    zgehrd.f  zstein.f \
claqge.f  dgtsv.f     dsptrd.f  slarnd.f    zgelq2.f  zstemr.f \
claqhb.f  dgtsvx.f    dsptrf.f  slarnv.f    zgelqf.f  zsteqr.f \
claqhe.f  dgttrf.f    dsptri.f  slarra.f    zgelsd.f  zsycon.f \
claqhp.f  dgttrs.f    dsptrs.f  slarrb.f    zgels.f   zsymv.f \
claqp2.f  dgtts2.f    dstebz.f  slarrc.f    zgelss.f  zsyr.f \
claqps.f  dhgeqz.f    dstedc.f  slarrd.f    zgelsx.f  zsyrfs.f \
claqr0.f  dhsein.f    dstegr.f  slarre.f    zgelsy.f  zsysv.f \
claqr1.f  dhseqr.f    dstein.f  slarrf.f    zgeql2.f  zsysvx.f \
claqr2.f  disnan.f    dstemr.f  slarrj.f    zgeqlf.f  zsytf2.f \
claqr3.f  dlabad.f    dsteqr.f  slarrk.f    zgeqp3.f  zsytrf.f \
claqr4.f  dlabrd.f    dsterf.f  slarrr.f    zgeqpf.f  zsytri.f \
claqr5.f  dlacn2.f    dstevd.f  slarrv.f    zgeqr2.f  zsytrs.f \
claqsb.f  dlacon.f    dstev.f   slartg.f    zgeqrf.f  ztbcon.f \
claqsp.f  dlacpy.f    dstevr.f  slartv.f    zgerfs.f  ztbrfs.f \
claqsy.f  dladiv.f    dstevx.f  slaruv.f    zgerq2.f  ztbtrs.f \
clar1v.f  dlae2.f     dsycon.f  slarzb.f    zgerqf.f  ztgevc.f \
clar2v.f  dlaebz.f    dsyevd.f  slarz.f     zgesc2.f  ztgex2.f \
clarcm.f  dlaed0.f    dsyev.f   slarzt.f    zgesdd.f  ztgexc.f \
clarfb.f  dlaed1.f    dsyevr.f  slas2.f     zgesvd.f  ztgsen.f \
clarf.f   dlaed2.f    dsyevx.f  slascl.f    zgesv.f   ztgsja.f \
clarfg.f  dlaed3.f    dsygs2.f  slasd0.f    zgesvx.f  ztgsna.f \
clarft.f  dlaed4.f    dsygst.f  slasd1.f    zgetc2.f  ztgsy2.f \
clarfx.f  dlaed5.f    dsygvd.f  slasd2.f    zgetf2.f  ztgsyl.f \
clargv.f  dlaed6.f    dsygv.f   slasd3.f    zgetrf.f  ztpcon.f \
clarnv.f  dlaed7.f    dsygvx.f  slasd4.f    zgetri.f  ztprfs.f \
clarrv.f  dlaed8.f    dsyrfs.f  slasd5.f    zgetrs.f  ztptri.f \
clartg.f  dlaed9.f    dsysv.f   slasd6.f    zggbak.f  ztptrs.f \
clartv.f  dlaeda.f    dsysvx.f  slasd7.f    zggbal.f  ztrcon.f \
clarzb.f  dlaein.f    dsytd2.f  slasd8.f    zgges.f   ztrevc.f \
clarz.f   dlaev2.f    dsytf2.f  slasda.f    zggesx.f  ztrexc.f \
clarzt.f  dlaexc.f    dsytrd.f  slasdq.f    zggev.f   ztrrfs.f \
clascl.f  dlag2.f     dsytrf.f  slasdt.f    zggevx.f  ztrsen.f \
claset.f  dlag2s.f    dsytri.f  slaset.f    zggglm.f  ztrsna.f \
clasr.f   dlags2.f    dsytrs.f  slasq1.f    zgghrd.f  ztrsyl.f \
classq.f  dlagtf.f    dtbcon.f  slasq2.f    zgglse.f  ztrti2.f \
claswp.f  dlagtm.f    dtbrfs.f  slasq3.f    zggqrf.f  ztrtri.f \
clasyf.f  dlagts.f    dtbtrs.f  slasq4.f    zggrqf.f  ztrtrs.f \
clatbs.f  dlagv2.f    dtgevc.f  slasq5.f    zggsvd.f  ztzrqf.f \
clatdf.f  dlahqr.f    dtgex2.f  slasq6.f    zggsvp.f  ztzrzf.f \
clatps.f  dlahr2.f    dtgexc.f  slasr.f     zgtcon.f  zung2l.f \
clatrd.f  dlahrd.f    dtgsen.f  slasrt.f    zgtrfs.f  zung2r.f \
clatrs.f  dlaic1.f    dtgsja.f  slassq.f    zgtsv.f   zungbr.f \
clatrz.f  dlaisnan.f  dtgsna.f  slasv2.f    zgtsvx.f  zunghr.f \
clatzm.f  dlaln2.f    dtgsy2.f  slaswp.f    zgttrf.f  zungl2.f \
clauu2.f  dlals0.f    dtgsyl.f  slasy2.f    zgttrs.f  zunglq.f \
clauum.f  dlalsa.f    dtpcon.f  slasyf.f    zgtts2.f  zungql.f \
cmach.f   dlalsd.f    dtprfs.f  slatbs.f    zhbevd.f  zungqr.f \
cpbcon.f  dlamch.f    dtptri.f  slatdf.f    zhbev.f   zungr2.f \
cpbequ.f  dlamrg.f    dtptrs.f  slatps.f    zhbevx.f  zungrq.f \
cpbrfs.f  dlaneg.f    dtrcon.f  slatrd.f    zhbgst.f  zungtr.f \
cpbstf.f  dlangb.f    dtrevc.f  slatrs.f    zhbgvd.f  zunm2l.f \
cpbsv.f   dlange.f    dtrexc.f  slatrz.f    zhbgv.f   zunm2r.f \
cpbsvx.f  dlangt.f    dtrrfs.f  slatzm.f    zhbgvx.f  zunmbr.f \
cpbtf2.f  dlanhs.f    dtrsen.f  slauu2.f    zhbtrd.f  zunmhr.f \
cpbtrf.f  dlansb.f    dtrsna.f  slauum.f    zhecon.f  zunml2.f \
cpbtrs.f  dlansp.f    dtrsyl.f  slazq3.f    zheevd.f  zunmlq.f \
cpocon.f  dlanst.f    dtrti2.f  slazq4.f    zheev.f   zunmql.f \
cpoequ.f  dlansy.f    dtrtri.f  sopgtr.f    zheevr.f  zunmqr.f \
cporfs.f  dlantb.f    dtrtrs.f  sopmtr.f    zheevx.f  zunmr2.f \
cposv.f   dlantp.f    dtzrqf.f  sorg2l.f    zhegs2.f  zunmr3.f \
cposvx.f  dlantr.f    dtzrzf.f  sorg2r.f    zhegst.f  zunmrq.f \
cpotf2.f  dlanv2.f    dzsum1.f  sorgbr.f    zhegvd.f  zunmrz.f \
cpotrf.f  dlapll.f    icmax1.f  sorghr.f    zhegv.f   zunmtr.f \
cpotri.f  dlapmt.f    ieeeck.f  sorgl2.f    zhegvx.f  zupgtr.f \
cpotrs.f  dlaptm.f    ilaenv.f  sorglq.f    zherfs.f  zupmtr.f \
cppcon.f  dlapy2.f    ilaver.f  sorgql.f    zhesv.f \
cppequ.f  dlapy3.f    iparmq.f  sorgqr.f    zhesvx.f \
cpprfs.f  dlaqgb.f    izmax1.f  sorgr2.f    zhetd2.f \
cppsv.f   dlaqge.f    lsame.f   sorgrq.f    zhetf2.f \
cppsvx.f  dlaqp2.f    lsamen.f  sorgtr.f    zhetrd.f \
clarnd.f  zlarnd.f


# Define name of library
libname = oomph_flapack

# Combine headers and sources
headers_and_sources =  $(sources)


# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = liboomph_flapack.la


# Sources that the library depends on:
#-------------------------------------
liboomph_flapack_la_SOURCES = $(headers_and_sources)
#Definitely no optimisation, but still warnings
liboomph_flapack_non_opt_la_FFLAGS="-O0"

# Automake doesn't have a way to build with no optimization
# See http://sources.redhat.com/ml/automake/2002-08/msg00082.html
# This is the automake-produced rule with $(FFLAGS) removed
slamch.lo: slamch.f
	$(LIBTOOL) --tag=F77 --mode=compile $(F77) $(FFLAGS_NO_OPT) -c -o slamch.lo `test -f 'slamch.f' || echo '$(srcdir)/'`slamch.f
dlamch.lo: dlamch.f
	$(LIBTOOL) --tag=F77 --mode=compile $(F77) $(FFLAGS_NO_OPT) -c -o dlamch.lo `test -f 'dlamch.f' || echo '$(srcdir)/'`dlamch.f


# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below  $(includedir)
include_HEADERS = $(libname).h

#Here's the subdirectory where the actual header files are placed
library_includedir=$(includedir)/oomph_flapack

#These are the header files that are to be placed in subdirectory: **NONE**
#library_include_HEADERS=$(headers)


# Combined header file **NONE**
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time
#$(libname).h: dummy_$(libname).h


# **NONE**
# When building the combined header, we also create a symbolic
# link to it from the directory above this one, so that
# it can be used by other libaries (before "make install" is run)
# in exactly the same way as after its installation in the include
# directory. (After the installation, the combined header lives
# in the directory above the actual individual header files)
#dummy_$(libname).h:  $(headers)
#	echo  $(libname) $(headers) > all_$(libname).aux
#	$(AWK) -f $(top_builddir)/bin/headers.awk < \
#	          all_$(libname).aux > $(libname).h
#	rm all_$(libname).aux
#	rm -f ../$(libname).h 
#	(cd .. && $(LN_S) $(libname)/$(libname).h $(libname).h)


#**NONE**
# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
#if SYMBOLIC_LINKS_FOR_HEADERS
#install-data-hook: 
#	(cd $(library_includedir) && rm -f $(headers) )
#	(echo "$(headers)" > include_files.list )
#	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
#	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
#	(rm -r include_files.list.aux include_files.list ) 
#else
install-data-hook: 
#endif

if ADD_LGFORTRAN_TO_OOMPH_LAPACK
liboomph_flapack_la_LIBADD = -lgfortran
endif