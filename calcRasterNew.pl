#!/usr/bin perl
#####################################################################
# This script calculates the size and angles at BPMA,BPMB and target
# in Hall A for a given raster powering.
# y. r. roblin december 2014.
# updated with new low disp optics December 2015.
####################################################################
# BPMA is first on diagnostic girder (IPM1H04A)
# BPMB is last on diagnostic girder (currently IPM1H04D)
use PDL;
use PDL::Matrix;
use Pod::Usage;

#
# quad powerings in Gauss directly from control system.
# (these default values are design values for 11023)
#($QKH01BDL,$QOH02BDL,$QMH02BDL,$QOH03BDL,$QOH03ABDL,$QAH04BDL)=(-58161,38184,-36785,38377,38660,-30801);
#($P,$QKH01BDL,$QOH02BDL,$QMH02BDL,$QOH03BDL,$QOH03ABDL,$QAH04BDL,$Ix,$Iy)=@ARGV;
($P,$Ix,$Iy)=@ARGV;

#($QKH01BDL,$QOH02BDL,$QMH02BDL,$QOH03BDL,$QOH03ABDL,$QAH04BDL)=(-58161,38184,-36785,38377,38660,-30801);

# LHRS 1722 RHRS 20843
($QKH01BDL,$QOH02BDL,$QMH02BDL,$QOH03BDL,$QOH03ABDL,$QAH04BDL)=(-2656.66, -5016.49, -3172.81, 3293.45, 3310.16, 3334.57);

#($QKH01BDL,$QOH02BDL,$QMH02BDL,$QOH03BDL,$QOH03ABDL,$QAH04BDL)=(-69438,39883,39797,-6260,32783,-56362);
my $ratio=$P/11023;
$QKH01BDL*=$ratio;
$QOH02BDL*=$ratio;
$QMH02BDL*=$ratio;
$QOH03BDL*=$ratio;
$QOH03ABDL*=$ratio;
$QAH04BDL*=$ratio;
#Ix and Iy are raster coil currents in Amps.
#
my $brho=$P*1e6/2.99792458e8;
#print "brho is $brho \n";

##################################################
# calculate quad transfer matrix
##################################################
sub quad {
    my ($BL,$L,$mref)=@_;
    my $G=$BL*1e2/$L/1e4;
    print "gradient is $G T/m \n";
    my $k=1.0*$G/$brho;
#    print "k1 is $k","\n";
    my $phi=$L/1e2*sqrt(abs($k));
#    print "phi is $phi \n";
    # for X defocusing quad
    do {
	set $mref, 0,0,cosh($phi);
	set $mref, 1,0,1.0/sqrt(abs($k))*sinh($phi);
	set $mref, 0,1,sqrt(abs($k))*sinh($phi);
	set $mref, 1,1,cosh($phi);
#
	set $mref, 2,2,cos($phi);
	set $mref, 3,2,1.0/sqrt(abs($k))*sin($phi);
	set $mref, 2,3,-sqrt(abs($k))*sin($phi);
	set $mref, 3,3,cos($phi);

    } if $k<0.0;
    # for X focusing quad
    do {
	set $mref, 0,0,cos($phi);
	set $mref, 1,0,1.0/sqrt($k)*sin($phi);
	set $mref, 0,1,-sqrt($k)*sin($phi);
	set $mref, 1,1,cos($phi);
#
	set $mref, 2,2,cosh($phi);
	set $mref, 3,2,1.0/sqrt(abs($k))*sinh($phi);
	set $mref, 2,3,sqrt(abs($k))*sinh($phi);
	set $mref, 3,3,cosh($phi);
    } if $k>0.0;
    
}
##################################################
# calculate drift transfer matrix
##################################################
sub drift {
    my ($L,$mref)=@_;
    set $mref, 0,0,1.0;
    set $mref, 1,0,$L/100.0;
    set $mref, 0,1,0;
    set $mref, 1,1,1.0;
#
    set $mref, 2,2,1.0;
    set $mref, 3,2,$L/100.0;
    set $mref, 2,3,0;
    set $mref, 3,3,1.0;
}

###################################################
#drifts
my $raster2H01=zeroes(4,4);
my $H01toOH02=zeroes(4,4);
my $OH02toMH02=zeroes(4,4);
my $MH02toOH03=zeroes(4,4);
my $OH03toOH03A=zeroes(4,4);
my $OH03AtoH04=zeroes(4,4);
my $H04toBPMA=zeroes(4,4);
my $BPMAtoBPMB=zeroes(4,4);
my $BPMBtoTgt=zeroes(4,4);
my $TgttoDiff=zeroes(4,4);
my $DifftoDump=zeroes(4,4);
#quads
my $QKH01=zeroes(4,4);
my $QOH02=zeroes(4,4);
my $QMH02=zeroes(4,4);
my $QOH03=zeroes(4,4);
my $QOH03A=zeroes(4,4);
my $QAH04=zeroes(4,4);
my $m2=zeroes(4,4);
my $drift=zeroes(4,4);
print "constructing beamline\n";
#
# construct beamline
# this is the layout as of december 17 2014.
drift(13210-12795,$raster2H01);
#drift(13210-13190,$raster2H01);
drift(13362-13240,$H01toOH02);
drift(13398-13422,$OH02toMH02);
drift(13495-13467,$MH02toOH03);
drift(13561-13531,$OH03toOH03A);
drift(14137-13597,$OH03AtoH04);
drift(14302-14167,$H04toBPMA);
#drift(14817-14302,$BPMAtoBPMB);
drift(514,$BPMAtoBPMB);
#drift(15005-14817,$BPMBtoTgt);
drift(123,$BPMBtoTgt);
drift(18255-15005,$TgttoDiff);
drift(20055-18255,$DifftoDump);
#TODO read these BDLs from the command line. 
# and also read the momentum
#
quad($QKH01BDL,30.0,$QKH01);
quad($QOH02BDL,36.22,$QOH02);
quad($QMH02BDL,45.05,$QMH02);
quad($QOH03BDL,36.22,$QOH03);
quad($QOH03ABDL,36.22,$QOH03A);
quad($QAH04BDL,30.0,$QAH04);
#my $ma=$raster2H01 x $QKH01;

my $RtoBPMA=zeroes(4,4);
$RtoBPMA=$H04toBPMA x $QAH04 x $OH03AtoH04 x $QOH03A x $OH03toOH03A x $QOH03 x $MH02toOH03 x $QMH02 x $OH02toMH02 x $QOH02 x $H01toOH02 x $QKH01 x $raster2H01;
my $RtoBPMB=zeroes(4,4);
$RtoBPMB=$BPMAtoBPMB x $RtoBPMA;
my $RtoTgt=zeroes(4,4);
$RtoTgt=$BPMBtoTgt x $RtoBPMB;
my $RtoDiff=zeroes(4,4);
$RtoDiff=$TgttoDiff x $RtoTgt;
my $RtoDump=zeroes(4,4);
$RtoDump=$DifftoDump x $RtoDiff;
#
#print $RtoBPMA," \n";
#print $RtoBPMB," \n";
#print $RtoTgt," \n";
#
# calculate the raster spot
#
#raster currents in Amp. 
# they can range from 0.0 to 40.0
#my $Ix=40.0;
#my $Iy=40.0;
# 86 G cm per Amp per coil (two coils per plane)
my $kickX=2*86/1e4*1e-2*$Ix/$brho;
my $kickY=2*86/1e4*1e-2*$Iy/$brho;
my $xspot=${RtoTgt}[1][0]*$kickX;
my $yspot=${RtoTgt}[3][2]*$kickY;
my $kickVec=pdl([0.0,$kickX,0.0,$kickY]);
my $sizesA=zeroes(1,4);
my $sizesB=zeroes(1,4);
my $sizesTgt=zeroes(1,4);
$sizesA=$RtoBPMA x $kickVec->transpose();
$sizesB=$RtoBPMB x $kickVec->transpose();
$sizesTgt=$RtoTgt x $kickVec->transpose();
$sizesDiff=$RtoDiff x $kickVec->transpose();
$sizesDump=$RtoDump x $kickVec->transpose();
print "size(x,x',y,y') at BPMA are $sizesA \n";
print "size(x,x',y,y') at BPMB are $sizesB \n";
print "size(x,x',y,y') at Tgt are $sizesTgt \n";
print "size(x,x',y,y') at Diff are $sizesDiff \n";
print "size(x,x',y,y') at Dump are $sizesDump \n";
