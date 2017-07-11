package kb_tophat2::kb_tophat2Client;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

kb_tophat2::kb_tophat2Client

=head1 DESCRIPTION


A KBase module: kb_tophat2


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => kb_tophat2::kb_tophat2Client::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 run_tophat2_app

  $returnVal = $obj->run_tophat2_app($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a kb_tophat2.TopHatInput
$returnVal is a kb_tophat2.TopHatResult
TopHatInput is a reference to a hash where the following keys are defined:
	input_ref has a value which is a kb_tophat2.obj_ref
	assembly_or_genome_ref has a value which is a kb_tophat2.obj_ref
	workspace_name has a value which is a string
	alignment_set_suffix has a value which is a string
	alignment_suffix has a value which is a string
	reads_condition has a value which is a string
	num_threads has a value which is an int
	read_mismatches has a value which is an int
	read_gap_length has a value which is an int
	read_edit_dist has a value which is an int
	min_intron_length has a value which is an int
	max_intron_length has a value which is an int
	min_anchor_length has a value which is an int
	report_secondary_alignments has a value which is a kb_tophat2.boolean
	no_coverage_search has a value which is a kb_tophat2.boolean
	library_type has a value which is a string
	preset_options has a value which is a string
obj_ref is a string
boolean is an int
TopHatResult is a reference to a hash where the following keys are defined:
	result_directory has a value which is a string
	reads_alignment_object_ref has a value which is a kb_tophat2.obj_ref
	report_name has a value which is a string
	report_ref has a value which is a string

</pre>

=end html

=begin text

$params is a kb_tophat2.TopHatInput
$returnVal is a kb_tophat2.TopHatResult
TopHatInput is a reference to a hash where the following keys are defined:
	input_ref has a value which is a kb_tophat2.obj_ref
	assembly_or_genome_ref has a value which is a kb_tophat2.obj_ref
	workspace_name has a value which is a string
	alignment_set_suffix has a value which is a string
	alignment_suffix has a value which is a string
	reads_condition has a value which is a string
	num_threads has a value which is an int
	read_mismatches has a value which is an int
	read_gap_length has a value which is an int
	read_edit_dist has a value which is an int
	min_intron_length has a value which is an int
	max_intron_length has a value which is an int
	min_anchor_length has a value which is an int
	report_secondary_alignments has a value which is a kb_tophat2.boolean
	no_coverage_search has a value which is a kb_tophat2.boolean
	library_type has a value which is a string
	preset_options has a value which is a string
obj_ref is a string
boolean is an int
TopHatResult is a reference to a hash where the following keys are defined:
	result_directory has a value which is a string
	reads_alignment_object_ref has a value which is a kb_tophat2.obj_ref
	report_name has a value which is a string
	report_ref has a value which is a string


=end text

=item Description

run_tophat2_app: run TopHat2 app

ref: https://ccb.jhu.edu/software/tophat/manual.shtml

=back

=cut

 sub run_tophat2_app
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_tophat2_app (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_tophat2_app:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_tophat2_app');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "kb_tophat2.run_tophat2_app",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_tophat2_app',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_tophat2_app",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_tophat2_app',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "kb_tophat2.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "kb_tophat2.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'run_tophat2_app',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method run_tophat2_app",
            status_line => $self->{client}->status_line,
            method_name => 'run_tophat2_app',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for kb_tophat2::kb_tophat2Client\n";
    }
    if ($sMajor == 0) {
        warn "kb_tophat2::kb_tophat2Client version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 boolean

=over 4



=item Description

A boolean - 0 for false, 1 for true.
@range (0, 1)


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 obj_ref

=over 4



=item Description

An X/Y/Z style reference


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 TopHatInput

=over 4



=item Description

required params:
input_ref: input reads object (Single/Paired_reads, reads_set, sample_set)
assembly_or_genome_ref: ref to Assembly, ContigSet, or Genome
workspace_name: the name of the workspace it gets saved to
alignment_set_suffix: suffix append to alignment set object name
alignment_suffix: suffix append to alignment object name

optional params:
reads_condition: condition associated with the input reads objec (ignored for sets of samples)
num_threads: number of processing threads
read_mismatches: read mismatch cutoff
read_gap_length: read gap cutoff
read_edit_dist: read edit cutoff
min_intron_length: minimum intron length
max_intron_length: maximum intron length
min_anchor_length: minimum anchor length
report_secondary_alignments: use this option to output secondary alignments
no_coverage_search: use this option to disable the coverage-based search for junctions
library_type: library type (fr-unstranded, fr-firststrand, fr-secondstrand)
preset_options: alignment preset options (b2-very-fast, b2-fast, b2-sensitive, b2-very-sensitive)

ref: https://ccb.jhu.edu/software/tophat/manual.shtml


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
input_ref has a value which is a kb_tophat2.obj_ref
assembly_or_genome_ref has a value which is a kb_tophat2.obj_ref
workspace_name has a value which is a string
alignment_set_suffix has a value which is a string
alignment_suffix has a value which is a string
reads_condition has a value which is a string
num_threads has a value which is an int
read_mismatches has a value which is an int
read_gap_length has a value which is an int
read_edit_dist has a value which is an int
min_intron_length has a value which is an int
max_intron_length has a value which is an int
min_anchor_length has a value which is an int
report_secondary_alignments has a value which is a kb_tophat2.boolean
no_coverage_search has a value which is a kb_tophat2.boolean
library_type has a value which is a string
preset_options has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
input_ref has a value which is a kb_tophat2.obj_ref
assembly_or_genome_ref has a value which is a kb_tophat2.obj_ref
workspace_name has a value which is a string
alignment_set_suffix has a value which is a string
alignment_suffix has a value which is a string
reads_condition has a value which is a string
num_threads has a value which is an int
read_mismatches has a value which is an int
read_gap_length has a value which is an int
read_edit_dist has a value which is an int
min_intron_length has a value which is an int
max_intron_length has a value which is an int
min_anchor_length has a value which is an int
report_secondary_alignments has a value which is a kb_tophat2.boolean
no_coverage_search has a value which is a kb_tophat2.boolean
library_type has a value which is a string
preset_options has a value which is a string


=end text

=back



=head2 TopHatResult

=over 4



=item Description

result_directory: folder path that holds all files generated by run_tophat2_app
reads_alignment_object_ref: generated Alignment/AlignmentSet object reference
report_name: report name generated by KBaseReport
report_ref: report reference generated by KBaseReport


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
result_directory has a value which is a string
reads_alignment_object_ref has a value which is a kb_tophat2.obj_ref
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
result_directory has a value which is a string
reads_alignment_object_ref has a value which is a kb_tophat2.obj_ref
report_name has a value which is a string
report_ref has a value which is a string


=end text

=back



=cut

package kb_tophat2::kb_tophat2Client::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
