package Tools::Taxonomy;

use strict;
use Bio::SeqIO;
use Data::Dumper;
use Logger::Logger;
use Cwd 'abs_path';
use base 'Exporter';
use Storable;


my $GG_DB = '/media/db/greengenes/gg_13_5_taxonomy.txt';


sub new {
  my ($class, %attrs) = @_;
  my $self = {};
  bless $self;
  $self->{_reduced_ranks} = { 'superkingdom' => 1,
                              'no rank' => 2,
                              'family' => 3,
                              'genus' => 4,
                              'species' => 5};

  $self->{_complete_ranks} = {'superkingdom' => 1,
                              'kingdom' => 2,
                              'subkingdom' => 3,
                              'superphylum' => 4,
                              'phylum' => 5,
                              'subphylum' => 6,
                              'superclass' => 7,
                              'class' => 8,
                              'subclass' => 9,
                              'superorder' => 10,
                              'order' => 11,
                              'suborder' => 12,
                              'superfamily' => 13,
                              'family' => 14,
                              'subfamily' => 15,
                              'tribe' => 16,
                              'subtribe' => 17,
                              'genus' => 18,
                              'subgenus' => 19,
                              'species' => 20,
                              'subspecies' => 21};

  if(defined($attrs{dbh})){
    $self->{dbh} = $attrs{dbh};
  }
  elsif(defined($attrs{gg})){
    $self->{gg} = 1;
    if($attrs{gg} != 1){
      _loadGgTaxonomy($self,$attrs{gg});
    }
    else{
      _loadGgTaxonomy($self,$GG_DB);
    }
  }
  else{
    $logger->warn('No database handle found.');
  }
  if(defined($attrs{_dbFormat})){
    if($attrs{_dbFormat} eq 'reduced'){
      $self->{_dbFormat} = 'reduced';
    }
  }

  return $self;
}

sub gi2taxonomy {
	my ($self,$gi,$algo)=@_;
	my $tax_id = retrieveTaxIdFromGi($self,$gi,$algo);
	my $taxo = _retrieveNCBITaxonomy($self,$tax_id);
	my $string = _getTaxonomyString($self,$taxo);
	return $string;
}

sub acc2taxonomy {
  my ($self,$acc,$algo)=@_;
  my $tax_id = retrieveTaxIdFromAcc($self,$acc,$algo);
  my $taxo = _retrieveNCBITaxonomy($self,$tax_id);
  my $string = _getTaxonomyString($self,$taxo);
	return $string;
}



sub _loadGgTaxonomy {
  my ($self,$ggTaxoFile) = @_;
  $logger->debug('Loading GreenGenes Taxonomy...');
  if(-e $ggTaxoFile . '.perl'){
    $logger->debug($ggTaxoFile . '.perl found, loading object...');
    $self->{gg_taxonomy} = _loadGgIndexFile($self,$ggTaxoFile . '.perl');
  }
  else{
    open(FILE,$ggTaxoFile) || $logger->logdie('GreenGenes taxonomy file not found. ' . $ggTaxoFile);
    my $c = { k => 'kingdom', p => 'phylum', c => 'class', o => 'order', f => 'family', g => 'genus', s => 'specie'};
    while(<FILE>){
      chomp;
      my @line = split(/\t/,$_);
      $line[1] =~ s/ //g;
      my @taxo = split(/;/,$line[1]);
      foreach my $s (@taxo){
        $s =~ /(\D)__(\S+)/;
        my $t = $1;
        my $u = $2;
        if(defined($c->{$t})){
          $self->{gg_taxonomy}->{$line[0]}->{$t} = $u;
        }
        else{
          $logger->error('Taxonomy type not defined.' . $t );
        }
      }
    }
    _writeGgObject($self,$self->{gg_taxonomy},$ggTaxoFile . '.perl' );
  }
  $logger->debug('Done.');
}


sub _loadGgIndexFile {
  my ($self,$file) = @_;
  my $index = retrieve($file);
  $logger->info('File '.$file." is now loaded\n");
  return $index;
}


sub _writeGgObject {

  my ($self,$index,$file) = @_;
  $logger->info('Writing index ('.(scalar keys %$index).' sequences) in file : '.$file."\n");
  store $index, $file;

  $logger->info('File '.$file." is now created\n");
  return $file;
}


sub retrieveTaxonomy {
  my ($self,$taxId) = @_;
  my $taxo;
  if(defined($self->{gg_taxonomy})){
    $taxo = _retrieveGreenGenesTaxonomy($self,$taxId);
  }
  else{
    $taxo = _retrieveNCBITaxonomy($self,$taxId);
  }
  return $taxo;
}


sub _retrieveGreenGenesTaxonomy {
  my ($self,$taxId) = @_;
  my $taxo;
  my $taxoMatch = {k => 'kingdom', p => 'phylum', c => 'class', o => 'order',f => 'family',g => 'genus',s =>'specie'};
  my @taxoOrder = ('k' , 'p' , 'c', 'p' ,'f' ,'g' ,'s' );
  $taxId =~ s/lcl\|//;
  if(defined($self->{gg_taxonomy}->{$taxId})){
    foreach my $t (@taxoOrder){
      if(defined($self->{gg_taxonomy}->{$taxId}->{$t})){
        push(@{$taxo},{ rank => $taxoMatch->{$t}, name => $self->{gg_taxonomy}->{$taxId}->{$t}, taxId => $taxId });
      }
      else{
        push(@{$taxo},{ rank => $taxoMatch->{$t}, name => 'unknown', taxId => $taxId });
      }
    }
  }
  else{
    $logger->error($taxId . ' not found in GreenGenes taxonomy.');
  }
  return $taxo;
}


sub _retrieveNCBITaxonomy {
  my ($self,$taxId) = @_;
  my $rank='';
  my $taxo;
  if(! defined($taxId)  || $taxId eq '' ){ return $taxo ; }
  if($taxId==1){
    unshift(@{$taxo}, {rank => 'root', name => 'root', 'taxId' => 1, 'parent_taxId' => 0 });
    return $taxo;
  }
  if($taxId<=0){
    unshift(@{$taxo}, {rank => 'unknown', name => 'unknown', 'taxId' => 0, 'parent_taxId' => 0 });
    return $taxo;
  }
  while($taxId != 1){
    if(! defined($taxId) || $taxId eq ''){
      if(defined($taxo)){
        return $taxo;
      }
      else{
        return 0;
      }
    }
    my $h;
    if(defined($self->{_storedTaxonomy}->{$taxId})){
      $logger->trace('already known tax_id');
      $h = $self->{_storedTaxonomy}->{$taxId};
      $logger->trace(Dumper $h);
    }
    else{
      my $cmd1 = 'SELECT * FROM nodes WHERE tax_id=' . $taxId;
      my $sth=$self->{dbh}->prepare($cmd1);
      $sth->execute();
      my @res = $sth->fetchrow_array();
      $rank = $res[2];

      my $parent_tax_id;
      if(defined($res[1])){
        $parent_tax_id = $res[1];
      }
      else{
        $parent_tax_id = undef;
      }
      $h = { 'rank' => $rank, 'name' => $self->retrieveNameFromTaxId($taxId), 'taxId' => $taxId, 'parent_taxId' => $parent_tax_id };
      $self->{_storedTaxonomy}->{$taxId} = $h;
    }
    if(defined($h->{parent_taxId})){
      $taxId=$h->{parent_taxId};
    }
    else{
      $logger->warn('Parent id not defined for tax_id: ' . $taxId);
      $h = {rank => 'unknown', name => 'unknown', 'taxId' => 0, 'parent_taxId' => 0 } ;
      $taxId=1;
    }

    unshift(@{$taxo}, $h);
  }
  if(defined $self->{_dbFormat} && $self->{_dbFormat} eq 'reduced'){
    my $reduced_taxo;
		foreach my $rank (sort {$self->{_reduced_ranks}->{$a} <=> $self->{_reduced_ranks}->{$b}} (keys(%{$self->{_reduced_ranks}}))){
			my $found=0;
			my $old_i=0;
			for(my $i=$old_i;$i<=$#{$taxo};$i++){
				if($rank eq 'no rank'){
					if($taxo->[$i-1]->{name} =~ /Viruses/ && $i==1){
						$found=1;
						push(@{$reduced_taxo},$taxo->[$i]);
						last;
					}
					$old_i=$i;
				}
				elsif($taxo->[$i]->{rank} eq $rank){
					push(@{$reduced_taxo},$taxo->[$i]);
					$found=1;
					$old_i=$i;
					last;
				}
			}
			if($found == 0){
				my $void = {rank => $rank, name => '', 'taxId' => undef, 'parent_taxId' => undef } ;
				push(@{$reduced_taxo},$void);
			}
		}
    $taxo = $reduced_taxo;
  }
  return $taxo;
}


sub retrieveTaxIdFromGi {
  my ($self,$gi,$algo) = @_;
  my $type;
  if(defined $algo && ($algo eq 'BLASTX' || $algo eq 'BLASTP' || $algo eq 'DIAMONDX' || $algo eq 'DIAMONDP')){
    $type = 'prot';
  }
  else{
    $type = 'nucl';
  }
  # my $NameQuery = 'SELECT tax_id FROM gi_taxid_' . $type . ' WHERE gi="' . $gi . '"';
  my $NameQuery = 'SELECT taxid FROM ' . $type . '_accession2taxid WHERE gi="' . $gi . '"';
  my $sth=$self->{dbh}->prepare($NameQuery);
  $sth->execute();
  my @res = $sth->fetchrow_array();
  if(! defined($res[0]) ){
    $logger->warn('tax_id not found for gi: ' . $gi . '.');
    return 0;
  }
  else{
    return $res[0];
  }
}


sub retrieveTaxIdFromAcc {
  my ($self,$acc,$algo) = @_;
  my $type;
  if ($acc =~ /(\w+)\.\d+/){
    $acc = $1
  }
  if(defined $algo && ($algo eq 'BLASTX' || $algo eq 'BLASTP' || $algo eq 'DIAMONDX' || $algo eq 'DIAMONDP')){
    $type = 'prot_accession2taxid';
  }
  else{
    $type = 'nucl_accession2taxid';
  }
  my $NameQuery = 'SELECT taxid FROM ' . $type . ' WHERE accession="' . $acc . '"';
  #my $NameQuery = 'SELECT taxid FROM ' . $type . ' WHERE accession="' . $acc . '"';
  $logger->debug($NameQuery);
  my $sth=$self->{dbh}->prepare($NameQuery);
  $sth->execute();
  my @res = $sth->fetchrow_array();
  if(! defined($res[0]) ){
    $logger->warn('tax_id not found for accession: ' . $acc . '.');
    return 0;
  }
  else{
    return $res[0];
  }
}


sub retrieveNameFromTaxId {
  my ($self,$taxId) = @_;
  my $nameQuery = 'SELECT name_txt FROM names WHERE tax_id="' . $taxId . '" AND name_class="scientific name"';
  my $sth=$self->{dbh}->prepare($nameQuery);
  $sth->execute();
  my @res = $sth->fetchrow_array();
  if(! defined($res[0])){
    return 0;
  }
  else{
    return $res[0];
  }
}


sub retrieveNameFromDescription {
  my ($self,$desc) = @_;
  if(!defined($desc)){
    return 'unknown';
  }
  if($desc =~ /^.*\[(.*)\]/){
    my $orga = $1;
    $logger->debug('Organism found in description : ' . $orga);
    return $orga;
  }
  else{
    $logger->warn('Organism not found in description string.');
    $logger->trace($desc);
    return 'unknown';
  }
}


sub retrieveTaxIdFromName {
  my ($self,$name) = @_;
  my $query = 'select tax_id from names where name_txt like "' . $name . '" limit 1';
  my $sth=$self->{dbh}->prepare($query);
  $sth->execute();
  my @res = $sth->fetchrow_array();
  if(! defined($res[0])){
    return 0;
  }
  else{
    $logger->debug('tax_id found from names :' . $res[0]);
    return $res[0];
  }
}


sub readCSVextended_regex {
  my ($self,$file,$separator,$regex) = @_ ;
  my $res = [];
  my @headers;

  open(CSV,$file) || $logger->logdie('Cannot open file ' . $file);
  while(<CSV>){
    chomp;
    if(/^#(.*)/){
      my $head = $1;
      $head =~ s/\"//g;
      @headers = split(/$separator/,$head);
    }
    else{
      if($_ !~ /$regex/){
        next;
      }
      my @line = split(/$separator/,$_);
      if(! @headers){
        @headers = (1 .. scalar @line);
      }
      my $m;
      for(my $i=0; $i<=$#line; $i++){
        $line[$i] =~ s/\"//g;
				# if(!defined($headers[$i])){
				# 	$logger->error(Dumper @headers);
				# 	$logger->error($i);
				# 	$logger->error(Dumper @line);
				# }
        $m->{$headers[$i]} = $line[$i];
      }
      push (@$res, $m);
    }
  }
  close CSV;
  $logger->info('ecsv file ' . $file . ' read, containing ' . scalar(@{$res}) . ' entries.');
  return wantarray ? ($res, \@headers) : $res;
}


sub readCSVextended {
  my ($self,$file,$separator) = @_ ;
  my $res = [];
  my @headers;
  open(CSV,$file) || $logger->logdie('Cannot open file ' . $file);
  while(<CSV>){
    chomp;
    if(/^#(.*)/){
      my $head = $1;
      $head =~ s/\"//g;
      @headers = split(/$separator/,$head);
    }
    else{
      my @line = split(/$separator/,$_);
      if(! @headers){
        @headers = (1 .. scalar @line);
      }
      my $m;
      for(my $i=0; $i<=$#line; $i++){
        $line[$i] =~ s/\"//g;
				# if(!defined($headers[$i])){
				# 	$logger->error(Dumper @headers);
				# 	$logger->error($i);
				# 	$logger->error(Dumper @line);
				# }
        $m->{$headers[$i]} = $line[$i];
      }
      push (@$res, $m);
    }
  }
  close CSV;
  $logger->info('ecsv file ' . $file . ' read, containing ' . scalar(@{$res}) . ' entries.');
  return wantarray ? ($res, \@headers) : $res;
}



sub _getNCBITaxonomy {
  my ($self,$m) = @_;
  if(! defined($m->{tax_id})){
    if( defined($m->{gi}) ){
      $m->{tax_id} = $self->retrieveTaxIdFromGi($m->{gi}, $m->{algo});
    }
    elsif(defined($m->{accession})){
      $m->{tax_id} = $self->retrieveTaxIdFromAcc($m->{accession}, $m->{algo});
    }
    else{
      $logger->warn('Neither gi nor accession defined.');
    }
  }

  if(! $m->{no_hit} && defined($m->{tax_id}) && $m->{tax_id} ne ''){
    if($m->{tax_id} == 0 && defined($m->{description}) && $m->{description} ne ''){
      $m->{organism} = $self->retrieveNameFromDescription( $m->{description} );
      if($m->{organism} ne 'unknown'){
        $m->{tax_id} = $self->retrieveTaxIdFromName($m->{organism});
        if($m->{tax_id} != 0){
          $m->{taxoTree} = $self->retrieveTaxonomy($m->{tax_id});
        }
        else{
          $m->{taxoTree} = [{rank => 'unknown', name => 'unknown', 'taxId' => 0 }];
        }
      }
      else{
        $m->{taxoTree} = [{rank => 'unknown', name => 'unknown', 'taxId' => 0 }];
      }
    }
    else{
      if($m->{tax_id}==1){
        $m->{taxoTree} = [{rank => 'root', name => 'root', taxId => 1 }];
      }
      elsif($m->{tax_id}<0){
        $m->{taxoTree} = [{rank => 'unknown', name => 'unknown', taxId => 0 }];
      }
      else{
        $m->{taxoTree} = $self->retrieveTaxonomy($m->{tax_id});
        $m->{organism} = $self->retrieveNameFromTaxId($m->{tax_id});
      }
    }
  }
  elsif($m->{no_hit}){
    $m->{taxoTree} = [{rank => 'unknown', name => 'No hit', 'taxId' => 0 }];
  }
  else{
    $m->{taxoTree} = [{rank => 'unknown', name => 'unknown', 'taxId' => 0 }];
  }
  $m->{taxonomy} = _getTaxonomyString($self,$m->{taxoTree});

  return $m->{taxoTree};
}


sub _getGreenGenesTaxonomy {
  my ($self,$m) = @_;
  $m->{tax_id} = $m->{gi} ;
  $m->{taxoTree} = $self->retrieveTaxonomy($m->{tax_id} );
  $m->{taxonomy} = _getTaxonomyString($self,$m->{taxoTree});
  $m->{organism} = _getGgOrganism($self,$m->{taxoTree});
}


sub _getGgOrganism {
  my ($self,$t) = @_;
  for(my $i=0;$i<$#{$t};$i++){
    if($t->[$#{$t}-$i]->{name} ne 'unknown'){
      return $t->[$#{$t}-$i]->{name};
    }
  }
}


sub _getTaxonomyString {
  my ($self,$t) = @_;
  my $s='';
  my $p=0;
  $s = join(';', map {$_->{'name'}} @$t) ;
  return $s;
}


sub buildFullTaxo {
  my ($self, $matches, $addQueryIdAsLeaf) = @_;
  my $tree ={};

  foreach my $x (@{$matches}){
    my $parent = $tree;

    my $taxoTree = _getNCBITaxonomy($self,$x);
    if($addQueryIdAsLeaf){
      push(@$taxoTree, {rank => 'unknown', name => $x->{query_id}, 'taxId' => 0 })
    }
    if(! $x->{nb_reads}){
      $x->{nb_reads}=1;
    }
    if(! $x->{nb_hit}){
      $x->{nb_hit}=1;
    }
    while(@$taxoTree){
      my $taxoRank = shift @$taxoTree;
      my $identity = $x->{percentIdentity};

      if(! defined $identity || $identity eq ''){
        $identity = 0;
      }
      if(defined($parent->{$taxoRank->{name}})){
        $parent->{$taxoRank->{name}}->{number} += $x->{nb_hit};
        $parent->{$taxoRank->{name}}->{identity}+= ($identity * $x->{nb_hit});
        $parent->{$taxoRank->{name}}->{number_r} += $x->{nb_reads};
        $parent->{$taxoRank->{name}}->{identity_r}+= ($identity * $x->{nb_reads});
      }
      else{
        $parent->{$taxoRank->{name}} = {
          number => $x->{nb_hit},
          number_r => $x->{nb_reads},
          rank   => $taxoRank->{rank},
          taxid  => $taxoRank->{taxId},
          identity => ($identity * $x->{nb_hit}),
          identity_r => ($identity * $x->{nb_reads})
        };
      }

      $parent = $parent->{$taxoRank->{name}};
    }
  }
  return $tree;
}

1;
