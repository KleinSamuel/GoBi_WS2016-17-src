package bamfiles;

/**
 * Used for mapper statistics as in assignment 3 task 2.
 */
public class Counter {

	private int NRP;
	private int MAPPED_NRP;
	private int MULTIMAPPED_NRP;
	private int TRANSCRIPTOMIC_RNP;
	private int MERGED_TR_NRP;
	private int INTRONIC_NRP;
	private int ANTISENSE_NRP;
	private int INTERGENIC_NRP;
	
	public Counter(){
		this.NRP = 0;
		this.MAPPED_NRP = 0;
		this.MULTIMAPPED_NRP = 0;
		this.TRANSCRIPTOMIC_RNP = 0;
		this.MERGED_TR_NRP = 0;
		this.INTRONIC_NRP = 0;
		this.ANTISENSE_NRP = 0;
		this.INTERGENIC_NRP = 0;
	}
	
	public String toString(){
		String s = "";
		s+= "Number Read Pairs:\t"+this.NRP+"\n";
		s+= "Mapped NRP:\t\t"+this.MAPPED_NRP+"\n";
		s+= "Multimapped NRP:\t"+this.MULTIMAPPED_NRP+"\n";
		s+= "Transcriptomic NRP:\t"+this.TRANSCRIPTOMIC_RNP+"\n";
		s+= "Merged_tr NRP:\t\t"+this.MERGED_TR_NRP+"\n";
		s+= "Intronic NRP:\t\t"+this.INTRONIC_NRP+"\n";
		s+= "Antisense NRP:\t\t"+this.ANTISENSE_NRP+"\n";
		s+= "Intergenic NRP:\t\t"+this.INTERGENIC_NRP+"\n";
		return s;
	}
	
	public int getNRP() {
		return NRP;
	}
	public int getMappedNRP() {
		return MAPPED_NRP;
	}
	public int getMultimappedNRP() {
		return MULTIMAPPED_NRP;
	}
	public int getTranscriptomicNRP() {
		return TRANSCRIPTOMIC_RNP;
	}
	public int getMergedTranscriptNRP() {
		return MERGED_TR_NRP;
	}
	public int getIntronicNRP() {
		return INTRONIC_NRP;
	}
	public int getAntisenseNRP() {
		return ANTISENSE_NRP;
	}
	public int getIntergenicNRP() {
		return INTERGENIC_NRP;
	}
	
	public void addNRP(){
		this.NRP++;
	}
	public void addMappedNRP(){
		this.MAPPED_NRP++;
	}
	public void addMultimappedNRP(){
		this.MULTIMAPPED_NRP++;
	}
	public void addTranscriptomicNRP(){
		this.TRANSCRIPTOMIC_RNP++;
	}
	public void addMergedTranscriptNRP(){
		this.MERGED_TR_NRP++;
	}
	public void addIntronicNRP(){
		this.INTRONIC_NRP++;
	}
	public void addAntisenseNRP(){
		this.ANTISENSE_NRP++;
	}
	public void addIntergenicNRP(){
		this.INTERGENIC_NRP++;
	}
	
}
