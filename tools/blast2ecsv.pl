#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use DBI;
use Tools::Blast;
use Tools::Taxonomy;
use Tools::Fasta;
use Pod::Usage;
use Logger::Logger;
use File::Basename;

my $VERSION = '1.2' ;
my $lastmodif = '2015-10-27' ;

my %blastOptions;
my $db = '/media/scratch/taxo/taxonomy.tmp.sqlite';
my $dbType = 'ncbi';
my $verbosity = 1;
my $parseMethod = 'global';
my $seqFile;
my $readNumber;
my @blastFiles;
my $outputFile = 'blast_ecsv.csv';
my $clean;
my $man;
my $help;
my $reduceTaxo;
my $vir_seq=0;


GetOptions(
  "b|blast=s"             => \@blastFiles,
  "o|output=s"            => \$outputFile,
  "pm|parsing_method=s"   => \$parseMethod,
  "t|type=s"				      => \$blastOptions{'algo'},
  "e|evalue=f"            => \$blastOptions{'max_evalue'},
  "s|score=i"             => \$blastOptions{'min_score'},
  "identity=s"            => \$blastOptions{'min_identity'},
  "qov|query_overlap=i"   => \$blastOptions{'min_query_overlap'},
  "hov|hit_overlap=i"     => \$blastOptions{'min_hit_overlap'},
  "if=s"                  => \$blastOptions{'format'},
  "dt|dbtype=s"           => \$dbType,
  "rm|remove_self_match"  => \$blastOptions{'remove_self_match'},
  "fhsp|first_hsp"        => \$blastOptions{'keep_first_hsp_only'},
  "fhit|first_hit"        => \$blastOptions{'keep_first_hit_only'},
  "pd"                    => \$blastOptions{'parse_description'},
  "rn|read_number=s"      => \$readNumber,
  "seq=s"                 => \$seqFile,
  "db=s"                  => \$db,
  "c|clean"               => \$clean,
  "r|reduce"              => \$reduceTaxo,
	"vs|vir_seq"            => \$vir_seq,
  "v|verbosity=i"         => \$verbosity,
  "h|help"                => \$help,
);
@blastFiles = split(/,/,join(',',@blastFiles));
if($help){
  &help;
}
if($verbosity > 1){
  Logger::Logger->changeMode($verbosity);
}


&main;


sub main {
  my $self={};
  bless $self;

  _set_options($self);
  my $fileInfos;
  my $fasta_tool;
  foreach my $file (@blastFiles){
    $self->{blast_tool} = Tools::Blast->new(file => $file, %blastOptions);

    my $hits = $self->deconvoluteBlastHits( $self->{blast_tool}->_readInputFile($self->{taxoTools}), $parseMethod );
    push(@$fileInfos, $hits);
  }
  # if($clean){
  # }
  $fileInfos = $self->mergeBlastCSV($fileInfos);

  my @blastExtentedHeaders = $self->getHeaders($blastOptions{'format'}, $readNumber);
  $self->printCSVExcel([$fileInfos], $outputFile, \@blastExtentedHeaders);
}


sub deconvoluteBlastHits {
  my ($self,$blastHits, $method) = @_;
  $logger->info('Methods ' . $method . '.');
  my $hits = [];
  if (!defined $method){
    $method ='local';
  }
  if($method =~ /^global$|^both$|^local$/){
    foreach my $match (@$blastHits){
      if($method eq 'global' || $method eq 'both'){
        push (@$hits, $match);
      }
      elsif($method eq 'local' || $method eq 'both'){
        if(defined($match->{no_hit}) && $method eq 'local'){
          push (@$hits, $match);
        }
        foreach my $hsp (@{$match->{'hsps'}}){
          push (@$hits, {%$match, %$hsp, 'nb_hsps' => 0});
        }
      }
    }
  }
  elsif($method eq 'lca'){
    @$hits = $self->_lca($blastHits);
  }
  return $hits;
}


sub _lca {
  my ($self,$blastHits)=@_;
  my $match_by_query_id={};
  for(my $i=0;$i<=$#{$blastHits};$i++){
    push(@{$match_by_query_id->{$blastHits->[$i]->{query_id}}->{index}},$i);
  }
  my @new_hit;
Q_ID:  foreach my $q_id (keys(%{$match_by_query_id})){
    foreach my $i (@{$match_by_query_id->{$q_id}->{index}}){
      if(! defined($blastHits->[$i]->{taxoTree})){
        next;
      }
      for(my $j=0;$j<=$#{$blastHits->[$i]->{taxoTree}};$j++){
        push(@{$match_by_query_id->{$q_id}->{ranks}->{ $blastHits->[$i]->{taxoTree}->[$j]->{rank} }->{names}}, $blastHits->[$i]->{taxoTree}->[$j]->{name});
      }
    }
    my $selected_name = '';
    if(defined($match_by_query_id->{$q_id}->{ranks})){
      $selected_name = $self->_computeRankHomogeneity( $match_by_query_id->{$q_id} );
      foreach my $i (@{$match_by_query_id->{$q_id}->{index}}){
        if(! defined($blastHits->[$i]->{taxonomy})){
          next;
        }
        if($blastHits->[$i]->{taxonomy} =~ /^(.*$selected_name);.*$/){
          push(@new_hit,{query_id => $q_id, taxonomy => $1, tax_id => $self->{taxoTools}->retrieveTaxIdFromName($selected_name)});
          next Q_ID;
        }
      }
    }
    else{
      push(@new_hit,{query_id => $q_id, taxonomy => 'no_hit', tax_id => 0});
    }
  }
  return @new_hit;
}


sub _computeRankHomogeneity {
  my ($self,$h)=@_;
  foreach my $rank (keys(%{$h->{ranks}})){
    foreach my $name (@{$h->{ranks}->{$rank}->{names}}){
      $h->{ranks}->{$rank}->{name_counter}->{$name}++;
    }
    $h->{ranks}->{$rank}->{total_name} = scalar(@{$h->{ranks}->{$rank}->{names}});
  }
  my $rank_order = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'];
  my $selected_name='';
  my $old_name='';
  foreach my $rank (@{$rank_order}){
    if(!defined($h->{ranks}->{$rank}->{total_name})){
      next;
    }
    if(($h->{ranks}->{$rank}->{total_name}/scalar(@{$h->{index}})) >= 0.6){
      foreach my $name (keys(%{$h->{ranks}->{$rank}->{name_counter}})){
        if(($h->{ranks}->{$rank}->{name_counter}->{$name}/$h->{ranks}->{$rank}->{total_name}) >= 0.8 ){
          $selected_name = $name;
        }
      }
    }
  }
  return $selected_name;
}


sub mergeBlastCSV {
  my ($self,$files) = @_;
  my $mergedFile= [];
	my $hash;
  for(my $i=0;$i<scalar(@{$files});$i++){
		for(my $j=0;$j<scalar(@{$files->[$i]});$j++){
			push(@{$hash->{$files->[$i]->[$j]->{query_id}}},$files->[$i]->[$j]);
		}
	}
	foreach my $query_id (keys(%{$hash})){
		if(defined($hash->{$query_id}->[0]->{no_hit})){
			my $found=0;
			for(my $i=1;$i<scalar(@{$hash->{$query_id}});$i++){
				if(defined($hash->{$query_id}->[$i]->{no_hit})){
					next;
				}
				else{
					push(@{$mergedFile},$hash->{$query_id}->[$i]);
					$found=1;
				}
			}
			if($found==0){
				push(@{$mergedFile},$hash->{$query_id}->[0]);
			}
		}
		else{
			push(@{$mergedFile},@{$hash->{$query_id}});
		}
	}
  return $mergedFile;
}


sub printCSVExcel {
  my ($self,$files, $outputFile, $headers, $separator) = @_;
  if( $headers && ! ref $headers ){
    $headers = [$headers];
  }
  if(!defined $separator){
    $separator = "\t";
  }
  open (OUTPUT, '>' . $outputFile) or $logger->logdie("Unable to open $outputFile");
  if($headers){
    print OUTPUT "#" . join($separator, @$headers)."\n";
  }

  foreach my $file ( @$files ){
    my $already_printed_sequences={};
    foreach my $line ( @$file ){
      if(! $headers || ! @$headers){
        $headers = [keys %$line];
        print OUTPUT "#" . join($separator, \@$headers)."\n";
      }
      my @fields;
      foreach my $field (@$headers){
        if($field eq 'sequence' && !defined($already_printed_sequences->{$line->{'query_id'}})){
          if($self->{_print_only_virus_seq} == 1){
            if($line->{'taxonomy'} =~ /Virus/ || $line->{'taxonomy'} =~ /Viroid/){
              #$logger->debug(ref($line->{'query_id'}));
              #if("$line->{'query_id'}" == "No"){next;}
              push(@fields, '"' . $self->{fasta_tool}->retrieveFastaSequence($line->{'query_id'})->{$line->{'query_id'}} . '"');
              $already_printed_sequences->{$line->{'query_id'}}=1;
            }
            else{
              push(@fields, '""');
              $already_printed_sequences->{$line->{'query_id'}}=1;
            }
          }
          else{
            if(length($self->{fasta_tool}->retrieveFastaSequence($line->{'query_id'})->{$line->{'query_id'}})>= 1){
              push(@fields, '"' . $self->{fasta_tool}->retrieveFastaSequence($line->{'query_id'})->{$line->{'query_id'}} . '"');
              $already_printed_sequences->{$line->{'query_id'}}=1;
            }
            else{
              push(@fields, '""');
              $already_printed_sequences->{$line->{'query_id'}}=1;
            }
          }
        }
        elsif($field eq 'nb_reads'){
          my $query_id = join( "", split("_d2b", $line->{'query_id'}) );
          if ($self->{read_numbers}->{$query_id}){
            push(@fields, '"' . $self->{read_numbers}->{$query_id} . '"');
          }else{
            push(@fields, '0');
          }
        }
        elsif($field eq 'query_length'){
          if(defined($line->{$field})){
            push(@fields, '"' . $line->{$field} . '"');
          }
          elsif(defined($seqFile)){
            push(@fields,'"' . length($self->{fasta_tool}->retrieveFastaSequence($line->{'query_id'})->{$line->{'query_id'}}) . '"');
          }
        }
        else{
          if (!defined $line->{$field}){
            $line->{$field} = "";
          }
          elsif(defined($line->{'no_hit'}) && $field !~ /algo|query_id|nb_reads|taxonomy/){
            $line->{$field} = "";
          }
          push(@fields, '"' . $line->{$field} . '"');
        }
      }
      print OUTPUT join($separator, @fields)."\n";
    }
  }
  close OUTPUT;
}


sub getHeaders {
  my ($self,$blastFormat,$readsPerContig) = @_;
  my @headers;
  if(defined $blastFormat){
    if($parseMethod eq 'lca'){
      @headers = ('query_id','tax_id','taxonomy');
    }
    else{
      if(defined($readsPerContig)){
        if($blastFormat eq 'xml'){
          @headers = ('algo', 'query_id', 'nb_reads', 'query_length', 'accession', 'description', 'organism', 'percentIdentity', 'nb_hsps', 'queryOverlap', 'hitOverlap', 'evalue', 'score', 'tax_id', 'taxonomy');
        }
        else{
          @headers = ('algo', 'query_id', 'nb_reads', 'accession', 'description', 'organism', 'percentIdentity', 'nb_hsps', 'evalue', 'score', 'tax_id', 'taxonomy');
        }
      }
      else{
        if($blastFormat eq 'xml'){
          @headers = ('algo', 'query_id', 'query_length', 'accession', 'description', 'organism', 'percentIdentity', 'evalue', 'tax_id', 'taxonomy');
        }
        else{
          @headers = ('algo', 'query_id', 'accession', 'description', 'organism', 'percentIdentity', 'nb_hsps', 'evalue', 'score', 'tax_id', 'taxonomy');
        }
      }
      if(defined($seqFile)){
        push(@headers,'sequence');
      }
    }
  }
  return @headers;
}


=head2 _readNumber

=head2

=head3 Description

Store reads number corresponding to each contig into a hash

=head3 Arguments

=over 4

=item

A TSV file containing the number of reads per contig (optional)

Format is :

Contig1 <TAB> ReadsNumer1

Contig2 <TAB> ReadsNumer2

...

=back

=head3 Returns

=over 4

=item

Hash reference containing contigs id as keys and number of reads as values

=back

=cut

sub _readNumber {
	my ($self,$readsFile) = @_;
	$logger->debug('Reading read number file...');
	my $reads;
	open(FILE, $readsFile) || $logger->logdie('File not found : ' . $readsFile);
	while(<FILE>){
		chomp;
		my @line = split(/\t/,$_);
		if(scalar(@line) != 2){
			$logger->logdie('Wrong file format. One line per entry, two columns per line.')
		}
		else{
			$reads->{$line[0]} = $line[1];
		}
	}
	$logger->debug('Done.');
	return $reads;
}


sub _set_options {
  my ($self) = @_;
  @{$self->{blastFiles}} = map{glob $_} @blastFiles;
  if(scalar(@{$self->{blastFiles}}) > 0){
    foreach my $file (@{$self->{blastFiles}}){
      if(! -e $file){
        $logger->error('Blast file not found. ' . $file);
        &help;
      }
    }
  }
  else{
    $logger->error('You must provide at least one Blast file.');
    &help;
  }


  my %taxonomyParams;
  if(! -e $db){
    $logger->error('Taxonomy sqlite database not found. ' . $db);
  }
  if($dbType eq 'ncbi'){
    $taxonomyParams{'dbh'} = DBI->connect("dbi:SQLite:dbname=$db","","");
  }
  elsif($dbType eq 'gg'){
    $taxonomyParams{'gg'} = 1;
  }
  else{
    $logger->error("Wrong database type. Must be ncbi or gg");
    &help();
  }
  if($reduceTaxo){
    $taxonomyParams{_dbFormat} = 'reduced';
  }
  else{
    $taxonomyParams{_dbFormat} = 'normal';
  }
  $self->{taxoTools} = Tools::Taxonomy->new(%taxonomyParams);

  if(defined($seqFile)){
    if(-e $seqFile){
      $self->{fasta_tool} = Tools::Fasta->new(file => $seqFile);
    }
    else{
      $logger->error('Fasta sequence file not found. ' . $seqFile),
    }
  }
  if(defined($readNumber)){
    if(-e $readNumber){
      $self->{read_numbers} = _readNumber($self,$readNumber);
    }
    else{
      $logger->error('Read Number file not found.' . $readNumber)
    }
  }
	if($vir_seq == 1){
		$logger->debug('vir_seq: will print only viral sequence.');
		$self->{_print_only_virus_seq} = 1;
	}
	else{
		$self->{_print_only_virus_seq} = 0;
	}

}


sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
### $prog $VERSION ###
#
# AUTHOR:     Sebastien THEIL
# VERSION:    $VERSION
# LAST MODIF: $lastmodif
# PURPOSE:    This script is used to parse blast+ output and produce a m8 output
#             extended with taxonomy information.
#

USAGE: perl $prog -i blast_1 -i blast_2 ... -i blast_n [OPTIONS]

          ### OPTIONS ###
          -b|blast                 <BLAST>  CSV or XML blast files
          -if                      <csv|xml>  Input  Blast file format (csv or xml). Default: csv
          -o|output                <ECSV>  Output ECSV file. Default: blast_ecsv.csv

          ### Blast parsing options ###
          -pm|parsing_method       <global|local|both>  Global: cumul hsps. Local: consider each hsp separetly. Both : keep cumuled hsp info and individual hsp infos. Default: local
          -fhsp|first_hsp          Treat only the first hsp
          -fhit|first_hit          Treat only the first hit
          -t|type                  <BLAST TYPE>  Force blast type detection (BLASTN|BLASTP|BLASTX|TBLASTX|TBLASTN|DIAMONDX|DIAMONDP). Usefull when using m8 blast format
          -e|evalue                <FLOAT>  E-value threshold
          -s|score                 <INTEGER>  Score threshold
          -identity                <FLOAT>  Identity threshold
          -qov                     <INTEGER>  Query overlap
          -hov                     <INTEGER>  Hit overlap
          -rm|remove_self_match    Remove self match
          -pd                      Parse description line to find query_id. Use it if blast has been launched without -parse_deflines.

          ### Additional inputs ###
          -rn|read_number          <READ NUMBER>  Tab-delimited file associating seqID with read number.
          -seq                     <FASTA>  Fasta sequence file. Usefull to retrieve sequence length when using m8 blast format.

          ### Taxonomy options ###
          -db                     <DATABASE>  Taxonomy database (Default:).
          -d|dbtype               <ncbi|gg>  Type of taxonomy database. ncbi = NCBI SQLite database, gg = greengene database; Default: ncbi.
          -r|reduce               Reduce taxonomy to 7 ranks : superkingdom;phylum;class;order;family;genus;species.

          ### Misc ###
          -c|clean                If multiple blast are provided and a query match, keep only one hit.
          -v|verbosity            Verbosity level
          -h|help                 Print this help and exit

EOF
exit(1) ;
}
