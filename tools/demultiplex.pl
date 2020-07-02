#!/usr/bin/perl

use strict;
use Data::Dumper;
use Logger::Logger;
use Getopt::Long;
use Tools::Fastq;
use Cwd 'abs_path', 'cwd' ;
use File::Basename ;
use String::Random;

my $VERSION = '3.2' ;
my $lastmodif = '2017-7-17' ;

my $pair1='';
my $pair2='';
my $indexFile='';
my $adapter='';
my $common={};
my $minLength=0;
my $quality=0;
my $pairsOutput = 'split';
my $force;
my $out='';
my $middle=0;
my $minIndexMatchPercent=0.4;
my $verbosity=1;
my $testMode=0;
my $index={};
my $polyA;
my $tmp_file_prefix='';
my $clean=0;



GetOptions(
            "-1=s"             => \$pair1,
            "-2=s"             => \$pair2,
            "indexFile=s"      => \$indexFile,
      "i|index=s%"       => \$index,
      "c|common=s%"      => \$common,
            "a|adapters=s"     => \$adapter,
            "l|len=i"          => \$minLength,
            "q|qual=i"         => \$quality,
            "po|pairsOutput=s" => \$pairsOutput,
            "f|force"          => \$force,
            "middle=i"         => \$middle,
            "o|outdir=s"       => \$out,
            "polyA"            => \$polyA,
            "test"             => \$testMode,
            "tmp_prefix=s"     => \$tmp_file_prefix,
            "clean"            => \$clean,
         "v|verbosity=i"    => \$verbosity,
);
sub help() {
  my $prog = basename($0) ;
  print STDERR <<EOF ;
  ### $prog $VERSION ###
  #
  # AUTHOR:     Sebastien Theil
  # VERSION:    $VERSION ($lastmodif)
  #
  # This script is used to demultiplex illumina reads using cutadapt
  # and trim off the adaptator sequence. It support paired end reads
  # and single reads.

  USAGE:
  Paired End Reads:
  $prog -1 R1.fq -2 R2.fq -i index.fa -c common.fa

  ### OPTIONS ###
  -1               Paired End reads 1.
  -2               Paired End reads 2.
  -indexFile       Index file sequences in fasta format.
  -i|index         Index tag and sequences in TAG=SEQUENCE format.
  -a|adapters      Illumina adapters sequences to trim.
 -c|common        Common sequence adapter in TAG=SEQUENCE format.
  -l|len           Minimum length of read to keep, after trimming.
  -q|qual          Minimum quality value.
  -pf|pairsFormat  Pairs format (split|merge). Separates pairs in two file or merge them.
  -po|pairsOutput  Pairs output (split|merge). Separates pairs depending on where they come. truePairs: same index. falsePairs: another index and singleton.
 -f|force
 -middle [1|2]    Search for common tag in the middle of the read. 1: trim the read. 2: exclude the read.
 -test            Test mode, just print the commands and don't execute them.
 -tmp_prefix      The prefix to name the tmp files. If none random one will be generated.
  -v|verbosity
  -h

EOF
  exit(1) ;
}

if($verbosity > 1){
    Logger::Logger->changeMode($verbosity);
}


&main;


sub main {
 my $self = {};
 bless $self;
 _setOptions($self);
 _launchPairedReadStep($self);
  _printStats($self);
 exit 0;
}


sub _step_01 {
  my ($self,$files,$tmp_file_prefix) = @_;
  _launchCutAdapt($self,$files->{1}, $self->{index}, $tmp_file_prefix . "_step.01_R1",'k','k','-g','0.1','1','0.8','1');
  _readInfoFile($self,$tmp_file_prefix . "_step.01_R1.info",'r1');
  _launchCutAdapt($self,$files->{2}, $self->{index}, $tmp_file_prefix . "_step.01_R2",'k','k','-g','0.1','1','0.8','1');
  _readInfoFile($self,$tmp_file_prefix . "_step.01_R2.info",'r2');
 return ($tmp_file_prefix . '_step.01_R1.out',$tmp_file_prefix . '_step.01_R2.out');
}


sub _step_02 {
  my ($self,$files,$tmp_file_prefix) = @_;
  _launchCutAdapt($self,$files->{1},$self->{_common},$tmp_file_prefix . "_step.02_R1",'k','k','-g','0.1',scalar(keys(%{$self->{_common}})),'0.7','2');
  _launchCutAdapt($self,$files->{2},$self->{_common},$tmp_file_prefix . "_step.02_R2",'k','k','-g','0.1',scalar(keys(%{$self->{_common}})),'0.7','2');
 if($self->{_clean} == 1 && $files->{1} ne $self->{readFiles}->{1}){
   push(@{$self->{_files_to_delete}},$files->{1});
   push(@{$self->{_files_to_delete}},$files->{2});
  }
 return ($tmp_file_prefix . '_step.02_R1.out',$tmp_file_prefix . '_step.02_R2.out');
}


sub _step_03 {
  my ($self,$files,$tmp_file_prefix) = @_;
  _launchCutAdapt($self,$files->{1},$self->{_common},$tmp_file_prefix . "_step.03_R1",'k','k','-g','0.2',scalar(keys(%{$self->{_common}})),'0.5','3');
  _launchCutAdapt($self,$files->{2},$self->{_common},$tmp_file_prefix . "_step.03_R2",'k','k','-g','0.2',scalar(keys(%{$self->{_common}})),'0.5','3');
 if($self->{_clean} == 1 && $files->{1} ne $self->{readFiles}->{1}){
   push(@{$self->{_files_to_delete}},$files->{1});
   push(@{$self->{_files_to_delete}},$files->{2});
  }
 return ($tmp_file_prefix . '_step.03_R1.out',$tmp_file_prefix . '_step.03_R2.out');
}


sub _step_04 {
  my ($self,$files,$tmp_file_prefix) = @_;
 my $h;
 foreach my $x (keys(%{$self->{illuminaAdapter}})){
  $h->{$x} = $self->{illuminaAdapter}->{$x};
  $h->{$x . '-REVCOMP'} = _reverseComplement($self,$self->{illuminaAdapter}->{$x});
 }
  $logger->info('Launch CutAdapt... step 04');
  _launchCutAdapt($self,$files->{1},$h,$tmp_file_prefix . "_step.04_R1",'k','k','-b','0.2','1','0.6','4');
  _launchCutAdapt($self,$files->{2},$h,$tmp_file_prefix . "_step.04_R2",'k','k','-b','0.2','1','0.6','4');
 if($self->{_clean} == 1 && $files->{1} ne $self->{readFiles}->{1}){
   push(@{$self->{_files_to_delete}},$files->{1});
   push(@{$self->{_files_to_delete}},$files->{2});
  }
 return ($tmp_file_prefix . '_step.04_R1.out',$tmp_file_prefix . '_step.04_R2.out');
}


sub _step_05 {
  my ($self,$files,$tmp_file_prefix) = @_;
 my $h;
 foreach my $x (keys(%{$self->{_common}})){
  $h->{$x} = $self->{_common}->{$x} ;
  $h->{$x . '-REVCOMP'} = _reverseComplement($self,$self->{_common}->{$x});
 }
 if($self->{middle_mid} == 1){
  _launchCutAdapt($self,$files->{1},$h,$tmp_file_prefix . "_step.05_R1",'k','k','-b','0.1','1','0.5','5');
   _launchCutAdapt($self,$files->{2},$h,$tmp_file_prefix . "_step.05_R2",'k','k','-b','0.1','1','0.5','5');
 }
  elsif($self->{middle_mid} == 2){
  _launchCutAdapt($self,$files->{1},$h,$tmp_file_prefix . "_step.05_R1",'k','d','-b','0.1','1','0.5','5');
   _launchCutAdapt($self,$files->{2},$h,$tmp_file_prefix . "_step.05_R2",'k','d','-b','0.1','1','0.5','5');
 }
 else{
  $logger->error('Wrong middle_mid value. ' . $self->{middle_mid});
  exit 1 ;
 }
 if($self->{_clean} == 1 && $files->{1} ne $self->{readFiles}->{1}){
   push(@{$self->{_files_to_delete}},$files->{1});
   push(@{$self->{_files_to_delete}},$files->{2});
  }
 return ($tmp_file_prefix . '_step.05_R1.out',$tmp_file_prefix . '_step.05_R2.out');
}


sub _step_06 {
  my ($self,$files,$tmp_file_prefix) = @_;
  my $h->{polyA} = 'AAAAAAAAAA';
 $h->{polyT} = 'TTTTTTTTTT';
  $logger->info('Launch CutAdapt... step 06');
  _launchCutAdapt($self,$files->{1}, $h, $tmp_file_prefix . "_step.06_R1",'k','k','-a','0','1','0.8','6');
  _launchCutAdapt($self,$files->{2}, $h, $tmp_file_prefix . "_step.06_R2",'k','k','-a','0','1','0.8','6');
  if($self->{_clean} == 1 && $files->{1} ne $self->{readFiles}->{1}){
   push(@{$self->{_files_to_delete}},$files->{1});
   push(@{$self->{_files_to_delete}},$files->{2});
  }
  return ($tmp_file_prefix . '_step.06_R1.out',$tmp_file_prefix . '_step.06_R2.out');
}


sub _reverseComplement {
 my ($self,$seq) = @_;
 my $revComp = reverse($seq);
 $revComp =~ tr/atcgnATCGN/tagcnTAGCN/;
 return $revComp;
}


sub _launchPairedReadStep {
 my ($self) = @_;
 my $i=0;
  my $files;
 my $last_files;
  $files->{1} = $self->{readFiles}->{1};
  $files->{2} = $self->{readFiles}->{2};

  # STEP 01: demultiplex from index file.
  if(defined($self->{index})){
    ($files->{1},$files->{2}) = _step_01($self,$files,$self->{tmp_file_prefix});

  }
  # STEP 02: 5' trimming common part.
 if(defined($self->{_common})) {
    ($files->{1},$files->{2}) = _step_02($self,$files,$self->{tmp_file_prefix});
  ($files->{1},$files->{2}) = _step_03($self,$files,$self->{tmp_file_prefix});
 }

  # STEP 04: 5' and 3' (reverse complement) trimming illumina adapters
 if(defined($self->{illuminaAdapter})) {
    ($files->{1},$files->{2}) = _step_04($self,$files,$self->{tmp_file_prefix});
 }

  # STEP 05: 3' trimming or excluding reads with common part (reverse complemented). (trimming chimeric reads)
 if(defined($self->{_common}) && $self->{middle_mid}>0) {
    ($files->{1},$files->{2}) = _step_05($self,$files,$self->{tmp_file_prefix});
 }

  if(defined($polyA)){
    ($files->{1},$files->{2}) = _step_06($self,$files,$self->{tmp_file_prefix});
  }

 my $last_R1_file = $files->{1};
 my $last_R2_file = $files->{2};
 $last_R1_file = abs_path($last_R1_file);
 $last_R2_file = abs_path($last_R2_file);

  _dispatchPairs($self,$last_R1_file,$last_R2_file);

  if($self->{_clean} == 1){
    foreach my $f (@{$self->{_files_to_delete}}){
      `rm $f`;
    }
  }
}


sub _get_pairs_and_singles {
  my ($self,$last_R1_file,$last_R2_file,$prefix)=@_;
  my @files = ($last_R1_file,$last_R2_file);
  $logger->info('Sorting pairs and singletons...');
  my $hash;
  my $old_hash;
  foreach my $f (@files){
    $logger->info('Reading fastq file: ' . $f);
    open(FILE,$f);
    while(<FILE>){
      if(/^@(\S+)\/([12])$/){
        if(defined($self->{dispatch})){
          if(defined($self->{dispatch}->{$prefix}->{$1})){
            $hash->{$1}->{$2}++;
          }
        }else{
          $hash->{$1}->{$2}++;
        }
      }
      if(/^>(\S+)\/([12])$/){
        if(defined($self->{dispatch})){
          if(defined($self->{dispatch}->{$prefix}->{$1})){
            $hash->{$1}->{$2}++;
          }
        }else{
          $hash->{$1}->{$2}++;
        }
      }
    }
    close FILE;
    $old_hash += $hash;
  }
  return $hash;
}

sub _create_sample_dir {
  my ($self,$name)=@_;
  my $path = $self->{_outDir} . '/' . $name;
  if(-e $path){
    $logger->warn('Directory ' . $path . ' already exists.');
    if(defined($self->{_force})){
      $logger->warn('Force mode...');
      `rm $path`;
    }
  }
  else{
    `mkdir $path`;
  }
  chdir $path;
  return $path;
}


sub _print_id_files {
  my ($self,$prefix,$hash)=@_;
  open(ID1,">$prefix" . "_r1.ids");
  open(ID2,">$prefix" . "_r2.ids");
  open(IDS,">$prefix" . "_s.ids");
  $logger->info("print id file");
  foreach my $r_id (sort(keys(%{$hash}))){
    my @keys = keys(%{$hash->{$r_id}});
  if(scalar(@keys) == 2){
   print ID1 $r_id . '/1' . "\n";
      print ID2 $r_id . '/2' . "\n";
  }
  else{
   foreach my $p (@keys){
    print IDS $r_id . '/' . $p . "\n";
   }
  }
  }
  close ID1;
  close ID2;
  close IDS;
}


sub _dispatchPairs {
  my ($self,$last_R1_file,$last_R2_file) = @_;

  $logger->info('Dispatching reads...');

  if(! defined($self->{index})){
    my $path = _create_sample_dir($self,$self->{tmp_file_prefix});
    my $hash = _get_pairs_and_singles($self,$last_R1_file,$last_R2_file,$self->{tmp_file_prefix});
    _print_id_files($self,$self->{tmp_file_prefix},$hash);

    _launch_subseq($self,$last_R1_file,$last_R2_file,$self->{tmp_file_prefix},$path);
  }
  else{
    foreach my $indexName (keys ( %{$self->{dispatch}}) ){
      $logger->info('Dispatching reads for index ' . $indexName);
      my $path = _create_sample_dir($self,$indexName);
      my $hash = _get_pairs_and_singles($self,$last_R1_file,$last_R2_file,$indexName);
      my $pwd = abs_path('.');
      $logger->info($pwd);

      _print_id_files($self,$indexName,$hash);

      _launch_subseq($self,$last_R1_file,$last_R2_file,$indexName,$path);
    }
  }
  chdir '..';
}


sub _launch_subseq {
  $logger->info("Launch subseq");
  my ($self,$last_R1_file,$last_R2_file,$prefix,$path)=@_;
  my $file1 = $path . '/' . $prefix . '_truePairs_r1.fq';
  my $file2 = $path . '/' . $prefix . '_truePairs_r2.fq';
  my $file_s = $path . '/' . $prefix . '_singletons.fq';

  my $seqtk_cmd = 'seqtk subseq ' . $last_R1_file . ' ' . $prefix . '_r1.ids > ' . $file1;
  $logger->debug($seqtk_cmd);
  $logger->info($seqtk_cmd);
  `$seqtk_cmd`;

  $seqtk_cmd = 'seqtk subseq ' . $last_R2_file . ' ' . $prefix . '_r2.ids > ' . $file2;
  $logger->debug($seqtk_cmd);
  $logger->info($seqtk_cmd);
  `$seqtk_cmd`;
  my $cmd = 'wc -l ' . $file1;
 $logger->debug($cmd);
  $logger->info($cmd);
  my ($nb_line_r1) = split(' ',`$cmd`);
 $logger->debug($nb_line_r1);
  $cmd = 'wc -l ' . $file2;
 $logger->debug($cmd);
  my ($nb_line_r2) = split(' ',`$cmd`);
 $logger->debug($nb_line_r2);
  if ($nb_line_r1 == $nb_line_r2){
    $self->{_final_stats}->{$prefix}->{pairs} = $nb_line_r1 / 4;
  }
  else{
    $logger->logdie('Read number between R1 and R2 are not he same. Something went wrong !!! r1:' . $nb_line_r1 . ' r2:' . $nb_line_r2);
  }
  $seqtk_cmd = 'seqtk subseq ' . $last_R1_file . ' ' . $prefix . '_s.ids > ' . $file_s;
  $logger->debug($seqtk_cmd);
  $logger->info($seqtk_cmd);
  `$seqtk_cmd`;
  $seqtk_cmd = 'seqtk subseq ' . $last_R2_file . ' ' . $prefix . '_s.ids >> ' . $file_s;
  $logger->debug($seqtk_cmd);
  $logger->info($seqtk_cmd);
  `$seqtk_cmd`;

  $cmd = 'wc -l ' . $file_s;
  my ($nb_reads_sing) = split(' ',`$cmd`);
  $self->{_final_stats}->{$prefix}->{singletons} = $nb_reads_sing / 4;
}


sub _printStats {
 my ($self) = @_;
  $logger->info('Printing statistics...');
 my $printed=0;
 my $stat_file = $self->{tmp_file_prefix} . '_demultiplex.stats.csv';
 open(STATS,">$stat_file") || $logger->logdie('Cannot create file ' . $stat_file);
  my @common_stats = ('total_reads','reads_ok','reads_with_adapters','reads_too_short','total_bases','bases_ok','bases_quality_trimmed');
  print STATS '#step,file,' . join(',',@common_stats) . "\n";
 foreach my $step (sort(keys(%{$self->{found}}))){
    foreach my $file (sort(keys(%{$self->{found}->{$step}}))){
      print STATS 'STEP' . $step . ',' . $file;
      foreach my $key (@common_stats){
        print STATS ',' . $self->{found}->{$step}->{$file}->{$key};
      }
      print STATS "\n";
    }
 }
  print STATS "\n";
  foreach my $step (sort(keys(%{$self->{found}}))){
    my $header_line = '#step,file';
    my $data_line = '';
    my %printed;
    foreach my $file (sort(keys(%{$self->{found}->{$step}}))){
      $data_line .= 'STEP' . $step . ',' . $file;
      foreach $index (sort(keys(%{$self->{found}->{$step}->{$file}->{index}}))){
        if(!defined($printed{$index})){
          $header_line .= ',' . $index . '(' . $self->{found}->{$step}->{$file}->{index}->{$index}->{seq} . ')';
        }
        $printed{$index}++;

        $data_line .= ',' . $self->{found}->{$step}->{$file}->{index}->{$index}->{nb};
      }
      $data_line .= "\n";
    }
    print STATS $header_line . "\n";
    print STATS $data_line . "\n";
  }
  print STATS '#index_id,pairs,singletons' . "\n";
  foreach $index (keys(%{$self->{_final_stats}})){
    print STATS $index . ',' . $self->{_final_stats}->{$index}->{pairs} . ',' . $self->{_final_stats}->{$index}->{singletons} . "\n";
  }
 close STATS;
}


sub _readInfoFile {
 my ($self,$file,$key) = @_;
 open(FILE,$file) || $logger->logdie('Info file not found ' . $file );
 $logger->info('Reading info file ' . $file);
 while(<FILE>){
  chomp;
  my @line = split(/\t/,$_);
  if(scalar(@line) <= 4){
   next;
  }
  my $match = { readName   => $line[0],
                errorNb    => $line[1],
                start      => $line[2],
                end        => $line[3],
                bLength    => length($line[4]),
                indexMatch => $line[5],
                aLength    => length($line[6]),
                indexName  => $line[7]
               };
  if($match->{readName} =~ /(\S+)\s\S+/){
   $match->{readName} = $1;
  }
    if($match->{readName} =~ /(\S+)\/([12])/){
      $match->{readName} = $1;
    }
  $self->{dispatch}->{ $match->{indexName} }->{ $match->{readName} } ++ ;
 }
}


sub _launchCutAdapt {
 my ($self,$file,$index,$prefix,$discardUntrimmed,$discardTrimmed,$indexLocation,$error,$count,$minIndexMatchPercent,$step_nb) = @_;
 if($step_nb == 1){
  if(-e $prefix . '.log' && -e $prefix . '.info'){
   $logger->info('Skipping cutadapt execution...');
   _readLogFile($self,$prefix . '.log',$step_nb);
   return 0;
  }
 }
 else{
  if(-e $prefix . '.log'){
   $logger->info('Skipping cutadapt execution...');
   _readLogFile($self,$prefix . '.log',$step_nb);
   return 0;
  }
 }

 my $max_length=999;
 $self->{_cutAdaptCmd} = 'cutadapt -f fastq' ;
 foreach my $p (keys(%{$index})){
  $self->{_cutAdaptCmd} .= ' ' . $indexLocation . ' ' . $p . '=' . $index->{$p} ;
  if(length($index->{$p}) < $max_length){
   $max_length = length($index->{$p});
  }
 }
 $self->{_cutAdaptCmd} .= ' -o ' . $prefix . '.out';
  if($step_nb == 1){
    $self->{_cutAdaptCmd} .= ' --info-file ' . $prefix . '.info';
  $self->{_cutAdaptCmd} .= ' --untrimmed-output=' . $prefix . '.no_mid.fastq';
  }
 if($discardUntrimmed eq 'd'){
  $self->{_cutAdaptCmd} .= ' --discard-untrimmed ' ;
 }
 if($discardTrimmed eq 'd'){
  $self->{_cutAdaptCmd} .= ' --discard-trimmed ' ;
 }
 if(defined($self->{_quality})){
  $self->{_cutAdaptCmd} .= ' -q ' . $self->{_quality} . ' --quality-base=' . $self->{_qualityBase} ;
 }
 if(defined($self->{_minLength})){
  $self->{_cutAdaptCmd} .= ' -m ' . $self->{_minLength} ;
 }
 if(defined($count) && $count > 0){
  $self->{_cutAdaptCmd} .= ' -n ' . $count;
 }
 my $overlapLength = int($max_length*$minIndexMatchPercent) ;
 $self->{_cutAdaptCmd} .= ' -O ' . $overlapLength;
 $self->{_cutAdaptCmd} .= ' -e ' . $error .
                          ' ' . $file;
 $self->{_cutAdaptCmd} .= ' > ' . $prefix . '.log';
 $logger->info($self->{_cutAdaptCmd});
 if(defined($self->{sge})){
  return $self->{_cutAdaptCmd};
 }
 else{
  if($self->{_testMode}==0){
   system($self->{_cutAdaptCmd});
   _readLogFile($self,$prefix . '.log',$step_nb);
  }
 }
}


sub _readLogFile {
 my ($self,$file,$step_nb) = @_;
 $logger->info('Reading log file ' . $file);
 open(LOG,$file) || $logger->logdie('Cannot open file ' . $file);
  my $seqName;
 while(<LOG>){
  chomp;
    if(/^Total reads processed:\s+(\S+)$/){
      $self->{found}->{$step_nb}->{$file}->{total_reads} = $1;
      $self->{found}->{$step_nb}->{$file}->{total_reads} =~ s/,//g;
    }
    if(/^Total basepairs processed:\s+(\S+)\sbp$/){
      $self->{found}->{$step_nb}->{$file}->{total_bases} = $1;
      $self->{found}->{$step_nb}->{$file}->{total_bases} =~ s/,//g;
    }
    if(/^Reads with adapters:\s+(\S+)\s\(\S+%\)$/){
      $self->{found}->{$step_nb}->{$file}->{reads_with_adapters} = $1;
      $self->{found}->{$step_nb}->{$file}->{reads_with_adapters} =~ s/,//g;
    }
    if(/^Reads written \(passing filters\):\s+(\S+)\s\(\S+%\)$/){
      $self->{found}->{$step_nb}->{$file}->{reads_ok} = $1;
      $self->{found}->{$step_nb}->{$file}->{reads_ok} =~ s/,//g;
    }
    if(/^Quality-trimmed:\s+(\S+)\sbp\s\(\S+%\)$/){
      $self->{found}->{$step_nb}->{$file}->{bases_quality_trimmed} = $1;
      $self->{found}->{$step_nb}->{$file}->{bases_quality_trimmed} =~ s/,//g;
    }
    if(/^Total written \(filtered\):\s+(\S+)\sbp\s\(\S+%\)$/){
      $self->{found}->{$step_nb}->{$file}->{bases_ok} = $1;
      $self->{found}->{$step_nb}->{$file}->{bases_ok} =~ s/,//g;
    }
    if(/^Reads that were too short:\s+(\S+)\s\(\S+\)$/){
      $self->{found}->{$step_nb}->{$file}->{reads_too_short} = $1;
      $self->{found}->{$step_nb}->{$file}->{reads_too_short} =~ s/,//g;
    }
    if(/^=== Adapter (\S+) ===$/){
      $seqName = $1;
    }
  if(/^Sequence: (\S+); Type: .*; Length: \d+; Trimmed: (\d+) times\./){
   $self->{found}->{$step_nb}->{$file}->{index}->{$seqName}->{nb} = $2;
   $self->{found}->{$step_nb}->{$file}->{index}->{$seqName}->{seq} = $1;
  }
 }
 close LOG;
}


sub _readIndexFile {
 my ($self,$file) = @_;
 open(FILE,$file) || $logger->logdie('File not found.');
 my $id ='';
 my $seq='';
 my $h;
 while(<FILE>){
  chomp;
  if(/^>(\S+)/){
   chomp;
   if($seq ne ''){
    $h->{$id} = $seq;
    $seq='';
   }
   $id=$1;
  }
  else{
   $seq.=$_;
  }
 }
 $h->{$id} = $seq;
 return $h;
}


#~ ------- TOOLS ------- ~#


sub _setOptions {
 my ($self) = @_;
 $self->{wd} = cwd;
 $self->{_qualityBase} = 33;

  if($clean==1){
    $self->{_clean}=1;
  }
  else{
    $self->{_clean}=0;
  }

 if($pair1 ne ''){
  if(-e $pair1){
   $self->{readFiles}->{1} = abs_path($pair1);
  }
  else{
   $logger->error('Pair1 file ' . $pair1 . ' do not exists.');
   &help;
  }
  $self->{_pairs} = 1;
 }

 if($pair2 ne ''){
  if(-e $pair2){
   $self->{readFiles}->{2} = abs_path($pair2);
  }
  else{
   $logger->error('Pair2 file ' . $pair2 . ' do not exists.');
   &help;
  }
  $self->{_pairs} = 1;
 }

 if($pair1 eq '' && $pair2 eq ''){
  $logger->error('You must provide at least one read file.');
  &help;
 }

 $self->{_quality} = $quality;
 $self->{_minLength} = $minLength;

 if($adapter ne ''){
  $self->{illuminaAdapter} = $self->_readIndexFile($adapter);
  $self->{_illuminaAdapterFile} = abs_path($indexFile);
 }

 if(scalar(keys(%{$common})) > 0){
  $self->{_common} = $common;
 }

 if(scalar(keys(%{$index})) > 0){
  $self->{index} = $index;
 }
 elsif($indexFile ne '' && scalar(keys(%{$index})) != 0){
  $self->{_indexFile} = abs_path($indexFile);
  $self->{index} = $self->_readIndexFile($self->{_indexFile});
 }
 else{
  $logger->warn('No index file provided.');
 }

 if($pairsOutput !~ /split|merge/ ){
  $logger->error('Wrong pairsOutput option. split|merge');
  &help;
 }
 else{
  $self->{_pairsFormat} = $pairsOutput;
 }

 if(defined($force)){
  $self->{_force} = 1;
 }

 if(defined($out)){
  $self->{_outDir} = abs_path($out);
 }
 else{
  $self->{_outDir} = abs_path('.');
 }

 if($middle != 0){
  $self->{middle_mid} = $middle;
 }
 else{
  $self->{middle_mid} = undef;
 }
 $self->{_minIndexMatchPercent} = $minIndexMatchPercent;
  if($testMode==1){
    $self->{_testMode} = 1 ;
  }
 else{
  $self->{_testMode} = 0 ;
 }
  if($tmp_file_prefix ne ''){
    $self->{tmp_file_prefix} = $tmp_file_prefix;
  }
  else{
    my $foo = new String::Random;
   $self->{tmp_file_prefix} = $foo->randregex('\w\w\w\w');
  }
}
