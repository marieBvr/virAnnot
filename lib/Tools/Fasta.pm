package Tools::Fasta;

use strict;
use warnings;
use Logger::Logger;
use Data::Dumper;
use Storable;

=head1 INDEXED FASTA RELATED METHODS

=head2

=head2 new

=head2

=head3 Description

Create a new Tools::Fasta object and index the FASTA file

=head3 Arguments

=over 4

=item

A hash of parameters.

Currently accepted keys are :

'file' => FASTA file path

=back

=head3 Returns

=over 4

=item

A Tools::Fasta object

=back

=cut

sub new {
	my ($class, %attrs) = @_;
	my $self = {};
	bless $self;
	if(defined($attrs{file})){
		$self->{file} = $attrs{file};
		open($self->{file_handle},$self->{file}) || $logger->logdie('Error opening file : '. $self->{file}.' : '.$!."\n");
		$self->indexFastaFile;
	}
	return $self;
}

=head2 indexFastaFile

=head2

=head3 Description

Index a FASTA file creating a hash reference with the following structure :

$index -> {seq_id} = {'id_begin_position' => integer, 'sequence_end_position' => integer, 'description' => string, 'position' => integer}

For each sequence id, the ">" symbol and all the text after space will be removed.

This cleaned id will be used as key for the index.

=head3 Arguments

=over 4

=item

None

=back

=head3 Returns

=over 4

=item

None

=back

=cut

sub indexFastaFile {
	my ($self) = @_;
  $logger->info('Indexing file : '.$self->{file}."\n");
  my $id;
  my $id_begin_position = 0;
  my $description;
	my $pos=0;
	my $fh = $self->{file_handle};
  while(my $line = <$fh>){
    if($line =~ /^>(\S+)\s?(.*)$/){
      $id = $1;
      $description = $2;
			$pos++;
      $self->{index}->{$id} = {'id_begin_position' => $id_begin_position, 'sequence_end_position' => 0, 'description' => $description, 'position' => $pos};
    }
    else{
      $self->{index}->{$id}{'sequence_end_position'} += length $line;
			# $logger->trace('TRACE: Indexing sequence' . $id . ' (id_begin_position : '. $self->{index}->{$id}{'id_begin_position'} . ', sequence_end_position : '. $self->{index}->{$id}{'sequence_end_position'} .') from ' . $self->{file} . "\n");
    }
    $id_begin_position = tell($fh);
  }
  $logger->info('File '.$self->{file}.' is now indexed (index contains '.(scalar keys %{$self->{index}})." sequences)\n");
}

=head2 loadFastaIndexFile

=head2

=head3 Description

Retrieve index from file using Storable module

=head3 Arguments

=over 4

An index file

=back

=head3 Returns

=over 4

=item

A hash reference corresponding to the index of the input FASTA file :

$index -> {seq_id} = {'id_begin_position' => integer, 'sequence_end_position' => integer, 'description' => string, 'position' => integer}

=back

=cut

sub loadFastaIndexFile {
	my ($self,$file) = @_;
  $self->{index} = retrieve($file);
	return $self->{index};
}

=head2 writeFastaIndexFile

=head2

=head3 Description

Write index to file using Storable module

=head3 Arguments

=over 4

=item

A hash reference corresponding to FASTA index.

=item

An output file path where to store the index.

=back

=head3 Returns

=over 4

=item

The output file path containing index

=back

=cut

sub writeFastaIndexFile {
	my ($self,$file) = @_;
  $logger->info('Writing index ('.(scalar keys %{$self->{index}}).' sequences) in file : '.$file."\n");
  store $self->{index}, $file;
  $logger->info('File '.$file." is now created\n");
	return $file;
}

=head2 retrieveFastaSequence

=head2

=head3 Description

Retrieve FASTA sequences using a list of ids

=head3 Arguments

=over 4

=item

A sequence id OR an array reference containing the list of sequences id to retrieve.

=back

=head3 Returns

=over 4

=item

A hash reference containing sequences id as keys and sequences as values

$data -> {seq_id} = sequence_corresponding_to_seq_id

=back

=cut

sub retrieveFastaSequence {
  my ($self,$ids) = @_;
  my $data;
  my $nbSequences = 0;
  my $fh = $self->{file_handle};
  if(! ref $ids){$ids = [$ids]}
  $logger->debug('Retrieving sequences of ids ['.join(', ', @$ids).'] from indexed file : '.$self->{file}."\n");
  foreach my $id (@$ids){
    my $cleanedId = $id;
    if($id =~ /^>?(\S*)/){$cleanedId = $1}
    #$logger->debug($cleanedId);
    $logger->debug('Retrieving informations of id ' . $cleanedId. " from index\n");
    if(exists $self->{index}->{$cleanedId}){
        # $logger->trace('DEBUG: id ' . $cleanedId . ' is present in index (id_begin_position : '. $self->{index}->{$cleanedId}{'id_begin_position'}. ', sequence_end_position : '. $self->{index}->{$cleanedId}{'sequence_end_position'}.")\n");
        seek($fh, $self->{index}->{$cleanedId}{'id_begin_position'}, 0);
        my $sequence = <$fh>;
        read($fh, $sequence, $self->{index}->{$cleanedId}{'sequence_end_position'});
        $sequence =~ s/\n//g;
      	$data->{$id} = $sequence;
        $nbSequences++;
    }
    else{
        # $logger->trace('DEBUG: id ' . $cleanedId. " not found in index\n");
    }
  }
  # $logger->trace($nbSequences.'/'.scalar(@$ids).' sequences has been retrieved from indexed file ' . $self->{file} . "\n");
  return $data;
}

=head2 retrieveFastaBlock

=head2

=head3 Description

Retrieve FASTA formatted sequences using a list of ids

=head3 Arguments

=over 4

=item

A sequence id OR an array reference containing the list of sequences id to retrieve.

=back

=head3 Returns

=over 4

=item

A scalar containing the sequences corresponding to ids in FASTA format

=back

=cut

sub retrieveFastaBlock {
	my ($self,$ids) = @_;
	my $data;
  my $nbSequences = 0;
	my $fh = $self->{file_handle};
  if(! ref $ids){$ids = [$ids]}
	$logger->debug('Retrieving fasta block of ' . scalar(@$ids) . ' ids from indexed file : '.$self->{file}."\n");
  foreach my $id (@$ids){
		my $cleanedId = $id;
		if($id =~ /^>?(\S*)/){$cleanedId = $1}
  	# $logger->trace('TRACE: retrieving informations of id ' . $cleanedId. " from index\n");
    if(exists $self->{index}->{$cleanedId}){
      # $logger->trace('TRACE: id ' . $cleanedId . ' is present in index (id_begin_position : '. $self->{index}->{$cleanedId}{'id_begin_position'}. ', sequence_end_position : '. $self->{index}->{$cleanedId}{'sequence_end_position'}.")\n");
      seek($fh, $self->{index}->{$cleanedId}{'id_begin_position'}, 0);
      my $blockId = <$fh>;
      read($fh, my $sequence, $self->{index}->{$cleanedId}{'sequence_end_position'});
    	$data .= $blockId . $sequence ;
      $nbSequences++;
      # $logger->trace('TRACE: fasta block of id ' . $cleanedId . ' is : ' . "\n" . $blockId.$sequence . "\n")
    }
    else{
			$logger->warn('WARN: id ' . $cleanedId. " not found in index\n");
		}
	}
  $logger->debug($nbSequences . '/' . scalar(@$ids) . ' fasta block has been retrieved from indexed file ' . $self->{file} . "\n");
	return $data;
}
1;
