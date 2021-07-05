#!usr/bin/env perl
use strict ;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Logger::Logger;
use Tools::Taxonomy;
use DBI;
use SQL::SplitStatement;


my $taxo_struct_dmp = 'taxonomyStructure.sql';
my $taxo_index_dmp = 'taxonomyIndex.sql';
my $data_acc_prot = 'prot.accession2taxid';
my $data_dead_acc_prot = 'dead_prot.accession2taxid';
my $data_dead_acc_nucl = 'dead_nucl.accession2taxid';
my $data_acc_wgs = 'nucl_wgs.accession2taxid';
my $data_acc_nucl = 'nucl.accession2taxid';
my $data_nodes = 'nodes.dmp';
my $data_names = 'names.dmp';
# my $gi_nucl = 'gi_taxid_nucl.dmp';
my $gi_prot = 'gi_taxid_prot.dmp';
my $dir = '.';
my $verbosity=1;


GetOptions(
  "acc_prot=s" => \$data_acc_prot,
  "acc_wgs=s"=> \$data_acc_wgs,
  "acc_nucl=s"=> \$data_acc_nucl,
  "dead_prot=s"=> \$data_dead_acc_prot,
  "dead_nucl=s"=> \$data_dead_acc_nucl,
  # "gi_nucl=s" => \$gi_nucl,
  "gi_prot=s" => \$gi_prot,
  "names=s"  => \$data_names,
  "nodes=s"  => \$data_nodes,
  "struct=s" => \$taxo_struct_dmp,
  "index=s"  => \$taxo_index_dmp,
  "dir=s"    => \$dir,
  "v=i"      => \$verbosity,
);

# `wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz ; gunzip taxdump.tar.gz; tar -xf taxdump.tar`;
# `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz ; gunzip prot.accession2taxid.gz`;
# `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz ; gunzip nucl_wgs.accession2taxid.gz`;
# `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz ; gunzip nucl_gb.accession2taxid.gz`;
# `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_prot.accession2taxid.gz ; gunzip dead_prot.accession2taxid.gz`;
# `wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz ; gunzip dead_wgs.accession2taxid.gz`;

Logger::Logger->changeMode($verbosity);


&main;


sub main {
  my $self={};
  bless $self;
  _set_options($self);

  $self->{_sep} = { names => '\t\|\t|\t\|$',
                    nodes => '\t\|\t|\t\|$',
                    'prot_accession2taxid' => '\t',
                    'nucl_accession2taxid' => '\t',
                    'gi_prot' => '\t'
                  };

  #
  my $db_path = $dir . '/taxonomy.tmp.sqlite';
  _create_sqlite_db($self,$db_path);
  my $dbh = DBI->connect("dbi:SQLite:dbname=$db_path","","");
  _insertingCSVDataInDatabase($self,$dbh,$self->{_data});
  $dbh->disconnect;
}

sub _insertingCSVDataInDatabase {
  my ($self,$dbh,$tablesDataFiles) = @_;
  $logger->info('Inserting tables into database...');
  foreach my $table (keys %{$tablesDataFiles}){
    $logger->info($table);
    my $sth = $dbh->column_info( undef, undef, $table, '%');
    my $ref = $sth->fetchall_arrayref;
    my @cols = map { $_->[3] } @$ref;
    my $req = "";
    $logger->debug("Inserting data in table $table ...\n");
    $dbh->{AutoCommit} = 0;
    if($table eq 'gi_nucl'){
      $req = "INSERT OR IGNORE INTO gi_taxid_nucl ( 'gi','tax_id' ) VALUES (?,?)"; 
    if($table eq 'gi_prot'){
      $req = "INSERT OR IGNORE INTO gi_taxid_prot ( 'gi','tax_id' ) VALUES (?,?)";  
    }else{
      $req = "INSERT OR IGNORE INTO $table ( ".join(',', map {"'".$_."'"} @cols)." ) VALUES (".join(',', map {'?'} @cols).")";  
    }
    
    $logger->info($req);
    # $sth = $dbh->prepare( "INSERT OR IGNORE INTO $table ( ".join(',', map {"'".$_."'"} @cols)." ) VALUES (".join(',', map {'?'} @cols).")" ) or $logger->logdie($dbh->errstr);
    $sth = $dbh->prepare( $req ) or $logger->logdie($dbh->errstr);

    my $separator;
    if(defined $self->{_sep}->{$table}){
      $separator = $self->{_sep}->{$table};
    }

    foreach my $file (@{$tablesDataFiles->{$table}}){
      $logger->info($file);
      open (DATA, $file) || $logger->logdie($file);
      while (<DATA>) {
        chomp;
        $sth->execute(grep {$_ !~ /^$separator$/} split (/($separator)/, $_)) or $logger->logdie($dbh->errstr);
      }
      close DATA;
    }

    $dbh->commit or $logger->logdie($dbh->errstr);
    $logger->debug("Insertion of data in table $table finished\n");
  }
}


sub _create_sqlite_db {
  my ($self,$file) = @_;
  $logger->info('Creating database.');
  if(-e $file){
    `mv $file $file.'_old'`;
  }
  if(! -e $file){
    `touch $file`;
    my $dbh = DBI->connect("dbi:SQLite:dbname=$file","","");
    _executeSQLFiles($self,$dbh,($self->{_taxo_struct_dmp}));
    $dbh->disconnect;
  }
  else{
    $logger->warn('Database already exists. Skip...')
  }
}


sub _executeSQLFiles {
  my ($self,$dbh,@sqlFiles) = @_;
  my $sql_splitter = SQL::SplitStatement->new;
  foreach my $file (@sqlFiles){
    $logger->debug('Reading sql file:' . $file);
    my $cmd;
    open (FILE, $file) or $logger->logdie("Unable to open the SQL file : $file\n");
    while( <FILE> ){
      $cmd.= $_;
    }
    close FILE;

    my @statements = $sql_splitter->split($cmd);
    foreach (@statements){
      $logger->debug('Executing sql cmd:');
      $logger->debug($_);
      $dbh-> do($_) or $logger->logdie($dbh->errstr);
    }
  }
}


sub _set_options {
  my ($self)=@_;
  if(-e $taxo_struct_dmp){
    $self->{_taxo_struct_dmp} = $taxo_struct_dmp;
  }
  else{
    $logger->error($taxo_struct_dmp . ' taxo_struct_dmp file not found.');
    &help;
  }


  if(-e $data_nodes){
    push(@{$self->{_data}->{nodes}},$data_nodes);
  }
  else{
    $logger->error($data_nodes . ' data_nodes file not found.');
    &help;
  }
  if(-e $data_names){
    push(@{$self->{_data}->{names}},$data_names);
  }
  else{
    $logger->error($data_names . ' data_names file not found.');
    &help;
  }
  if(-e $data_acc_prot){
    push(@{$self->{_data}->{prot_accession2taxid}},$data_acc_prot);
  }
  else{
    $logger->error($data_acc_prot . ' data_acc_prot file not found.');
    &help;
  }
  if(-e $data_dead_acc_prot){
    push(@{$self->{_data}->{prot_accession2taxid}},$data_dead_acc_prot);
  }
  else{
    $logger->error($data_dead_acc_prot . ' data_acc_prot file not found.');
    &help;
  }
  if(-e $data_dead_acc_nucl){
    push(@{$self->{_data}->{nucl_accession2taxid}},$data_dead_acc_nucl);
  }
  else{
    $logger->error($data_dead_acc_nucl . ' data_dead_acc_nucl file not found.');
    &help;
  }
  if(-e $data_acc_nucl){
    push(@{$self->{_data}->{nucl_accession2taxid}},$data_acc_nucl);
  }
  else{
    $logger->error($data_acc_nucl . ' data_acc_nucl file not found.');
    &help;
  }
  if(-e $data_acc_wgs){
    push(@{$self->{_data}->{nucl_accession2taxid}},$data_acc_wgs);
  }
  else{
    $logger->error($data_acc_wgs . ' data_acc_wgs file not found.');
    &help;
  }
  # if(-e $gi_nucl){
  #   push(@{$self->{_data}->{gi_nucl}},$gi_nucl);
  # }
  # else{
  #   $logger->error($gi_nucl . ' gi_nucl file not found.');
  #   &help;
  # }
  if(-e $gi_prot){
    push(@{$self->{_data}->{gi_prot}},$gi_prot);
  }
  else{
    $logger->error($gi_prot . ' gi_nucl file not found.');
    &help;
  }
  # if(-e $data_acc_gss){
  #   push(@{$self->{_data}->{nucl_accession2taxid}},$data_acc_gss);
  # }
  # else{
  #   $logger->error($data_acc_gss . ' data_acc_wgs file not found.');
  #   &help;
  # }
  # if(-e $data_acc_gb){
  #   push(@{$self->{_data}->{nucl_accession2taxid}},$data_acc_gb);
  # }
  # else{
  #   $logger->error($data_acc_gb . ' data_acc_gb file not found.');
  #   &help;
  # }
}


  sub help {
    my $prog = basename($0);
    print STDERR <<EOF ;
    #### $prog ####
    #
    # AUTHOR:     Sebastien THEIL and Marie LEFEBVRE
    # LAST MODIF: 07/02/2020
    # PURPOSE:    This script is used to load NCBI taxonomy file into a SQLite database.

    USAGE:
          $prog -struct taxonomyStructure.sql -index taxonomyIndex.sql -acc_prot acc2taxid.prot -acc_nucl acc2taxid.nucl -names names.dmp -nodes nodes.dmp -gi_nucl gi_taxid_nucl.dmp -gi_prot gi_taxid_prot.dmp

          ### OPTIONS ###
          -struct     <path>   taxonomyStructure.sql path. (Default: $taxo_struct_dmp)
          -index      <path>   taxonomyIndex.sql path. (Default: $taxo_index_dmp)
          -acc_prot   <path>   prot.accession2taxid. (Default: $data_acc_prot)
          -acc_nucl   <path>   nucl.accession2taxid. (Default: $data_acc_wgs)
          -names      <path>   names.dmp file. (Default: $data_names)
          -nodes      <path>   nodes.dmp file. (Default: $data_nodes)
          -gi_prot    <path>   gi_taxid_prot.dmp file (Default: $gi_prot)
          -v          <int>      Verbosity level. (0 -> 4).
    EOF
    exit(1);
  }
}
