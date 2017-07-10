
package us.kbase.kbtophat2;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: TopHatInput</p>
 * <pre>
 * required params:
 * input_ref: input reads object (Single/Paired_reads, reads_set, sample_set)
 * assembly_or_genome_ref: ref to Assembly, ContigSet, or Genome
 * workspace_name: the name of the workspace it gets saved to
 * alignment_object_name: output Alignment or AlignmentSet object name
 * optional params:
 * reads_condition: condition associated with the input reads objec (ignored for sets of samples)
 * num_threads: number of processing threads
 * read_mismatches: read mismatch cutoff
 * read_gap_length: read gap cutoff
 * read_edit_dist: read edit cutoff
 * min_intron_length: minimum intron length
 * max_intron_length: maximum intron length
 * min_anchor_length: minimum anchor length
 * report_secondary_alignments: use this option to output secondary alignments
 * no_coverage_search: use this option to disable the coverage-based search for junctions
 * library_type: library type (fr-unstranded, fr-firststrand, fr-secondstrand)
 * preset_options: alignment preset options (b2-very-fast, b2-fast, b2-sensitive, b2-very-sensitive)
 * ref: https://ccb.jhu.edu/software/tophat/manual.shtml
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "input_ref",
    "assembly_or_genome_ref",
    "workspace_name",
    "alignment_object_name",
    "reads_condition",
    "num_threads",
    "read_mismatches",
    "read_gap_length",
    "read_edit_dist",
    "min_intron_length",
    "max_intron_length",
    "min_anchor_length",
    "report_secondary_alignments",
    "no_coverage_search",
    "library_type",
    "preset_options"
})
public class TopHatInput {

    @JsonProperty("input_ref")
    private String inputRef;
    @JsonProperty("assembly_or_genome_ref")
    private String assemblyOrGenomeRef;
    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("alignment_object_name")
    private String alignmentObjectName;
    @JsonProperty("reads_condition")
    private String readsCondition;
    @JsonProperty("num_threads")
    private Long numThreads;
    @JsonProperty("read_mismatches")
    private Long readMismatches;
    @JsonProperty("read_gap_length")
    private Long readGapLength;
    @JsonProperty("read_edit_dist")
    private Long readEditDist;
    @JsonProperty("min_intron_length")
    private Long minIntronLength;
    @JsonProperty("max_intron_length")
    private Long maxIntronLength;
    @JsonProperty("min_anchor_length")
    private Long minAnchorLength;
    @JsonProperty("report_secondary_alignments")
    private Long reportSecondaryAlignments;
    @JsonProperty("no_coverage_search")
    private Long noCoverageSearch;
    @JsonProperty("library_type")
    private String libraryType;
    @JsonProperty("preset_options")
    private String presetOptions;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("input_ref")
    public String getInputRef() {
        return inputRef;
    }

    @JsonProperty("input_ref")
    public void setInputRef(String inputRef) {
        this.inputRef = inputRef;
    }

    public TopHatInput withInputRef(String inputRef) {
        this.inputRef = inputRef;
        return this;
    }

    @JsonProperty("assembly_or_genome_ref")
    public String getAssemblyOrGenomeRef() {
        return assemblyOrGenomeRef;
    }

    @JsonProperty("assembly_or_genome_ref")
    public void setAssemblyOrGenomeRef(String assemblyOrGenomeRef) {
        this.assemblyOrGenomeRef = assemblyOrGenomeRef;
    }

    public TopHatInput withAssemblyOrGenomeRef(String assemblyOrGenomeRef) {
        this.assemblyOrGenomeRef = assemblyOrGenomeRef;
        return this;
    }

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public TopHatInput withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("alignment_object_name")
    public String getAlignmentObjectName() {
        return alignmentObjectName;
    }

    @JsonProperty("alignment_object_name")
    public void setAlignmentObjectName(String alignmentObjectName) {
        this.alignmentObjectName = alignmentObjectName;
    }

    public TopHatInput withAlignmentObjectName(String alignmentObjectName) {
        this.alignmentObjectName = alignmentObjectName;
        return this;
    }

    @JsonProperty("reads_condition")
    public String getReadsCondition() {
        return readsCondition;
    }

    @JsonProperty("reads_condition")
    public void setReadsCondition(String readsCondition) {
        this.readsCondition = readsCondition;
    }

    public TopHatInput withReadsCondition(String readsCondition) {
        this.readsCondition = readsCondition;
        return this;
    }

    @JsonProperty("num_threads")
    public Long getNumThreads() {
        return numThreads;
    }

    @JsonProperty("num_threads")
    public void setNumThreads(Long numThreads) {
        this.numThreads = numThreads;
    }

    public TopHatInput withNumThreads(Long numThreads) {
        this.numThreads = numThreads;
        return this;
    }

    @JsonProperty("read_mismatches")
    public Long getReadMismatches() {
        return readMismatches;
    }

    @JsonProperty("read_mismatches")
    public void setReadMismatches(Long readMismatches) {
        this.readMismatches = readMismatches;
    }

    public TopHatInput withReadMismatches(Long readMismatches) {
        this.readMismatches = readMismatches;
        return this;
    }

    @JsonProperty("read_gap_length")
    public Long getReadGapLength() {
        return readGapLength;
    }

    @JsonProperty("read_gap_length")
    public void setReadGapLength(Long readGapLength) {
        this.readGapLength = readGapLength;
    }

    public TopHatInput withReadGapLength(Long readGapLength) {
        this.readGapLength = readGapLength;
        return this;
    }

    @JsonProperty("read_edit_dist")
    public Long getReadEditDist() {
        return readEditDist;
    }

    @JsonProperty("read_edit_dist")
    public void setReadEditDist(Long readEditDist) {
        this.readEditDist = readEditDist;
    }

    public TopHatInput withReadEditDist(Long readEditDist) {
        this.readEditDist = readEditDist;
        return this;
    }

    @JsonProperty("min_intron_length")
    public Long getMinIntronLength() {
        return minIntronLength;
    }

    @JsonProperty("min_intron_length")
    public void setMinIntronLength(Long minIntronLength) {
        this.minIntronLength = minIntronLength;
    }

    public TopHatInput withMinIntronLength(Long minIntronLength) {
        this.minIntronLength = minIntronLength;
        return this;
    }

    @JsonProperty("max_intron_length")
    public Long getMaxIntronLength() {
        return maxIntronLength;
    }

    @JsonProperty("max_intron_length")
    public void setMaxIntronLength(Long maxIntronLength) {
        this.maxIntronLength = maxIntronLength;
    }

    public TopHatInput withMaxIntronLength(Long maxIntronLength) {
        this.maxIntronLength = maxIntronLength;
        return this;
    }

    @JsonProperty("min_anchor_length")
    public Long getMinAnchorLength() {
        return minAnchorLength;
    }

    @JsonProperty("min_anchor_length")
    public void setMinAnchorLength(Long minAnchorLength) {
        this.minAnchorLength = minAnchorLength;
    }

    public TopHatInput withMinAnchorLength(Long minAnchorLength) {
        this.minAnchorLength = minAnchorLength;
        return this;
    }

    @JsonProperty("report_secondary_alignments")
    public Long getReportSecondaryAlignments() {
        return reportSecondaryAlignments;
    }

    @JsonProperty("report_secondary_alignments")
    public void setReportSecondaryAlignments(Long reportSecondaryAlignments) {
        this.reportSecondaryAlignments = reportSecondaryAlignments;
    }

    public TopHatInput withReportSecondaryAlignments(Long reportSecondaryAlignments) {
        this.reportSecondaryAlignments = reportSecondaryAlignments;
        return this;
    }

    @JsonProperty("no_coverage_search")
    public Long getNoCoverageSearch() {
        return noCoverageSearch;
    }

    @JsonProperty("no_coverage_search")
    public void setNoCoverageSearch(Long noCoverageSearch) {
        this.noCoverageSearch = noCoverageSearch;
    }

    public TopHatInput withNoCoverageSearch(Long noCoverageSearch) {
        this.noCoverageSearch = noCoverageSearch;
        return this;
    }

    @JsonProperty("library_type")
    public String getLibraryType() {
        return libraryType;
    }

    @JsonProperty("library_type")
    public void setLibraryType(String libraryType) {
        this.libraryType = libraryType;
    }

    public TopHatInput withLibraryType(String libraryType) {
        this.libraryType = libraryType;
        return this;
    }

    @JsonProperty("preset_options")
    public String getPresetOptions() {
        return presetOptions;
    }

    @JsonProperty("preset_options")
    public void setPresetOptions(String presetOptions) {
        this.presetOptions = presetOptions;
    }

    public TopHatInput withPresetOptions(String presetOptions) {
        this.presetOptions = presetOptions;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((((((((((((((((((((((("TopHatInput"+" [inputRef=")+ inputRef)+", assemblyOrGenomeRef=")+ assemblyOrGenomeRef)+", workspaceName=")+ workspaceName)+", alignmentObjectName=")+ alignmentObjectName)+", readsCondition=")+ readsCondition)+", numThreads=")+ numThreads)+", readMismatches=")+ readMismatches)+", readGapLength=")+ readGapLength)+", readEditDist=")+ readEditDist)+", minIntronLength=")+ minIntronLength)+", maxIntronLength=")+ maxIntronLength)+", minAnchorLength=")+ minAnchorLength)+", reportSecondaryAlignments=")+ reportSecondaryAlignments)+", noCoverageSearch=")+ noCoverageSearch)+", libraryType=")+ libraryType)+", presetOptions=")+ presetOptions)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
