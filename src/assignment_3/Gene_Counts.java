package assignment_3;

public class Gene_Counts {
	
	private int NRPwithinGenReg;
	private int NRPwithinGenRegWitoutOtherGenReg;
	private int NRPintronic;
	private int NRPtranscriptomic;
	private int NRPmergedTr;
	private int NRPforHighestTranscript;
	
	
	public Gene_Counts(){
		this.NRPwithinGenReg = 0;
		this.NRPwithinGenRegWitoutOtherGenReg = 0;
		this.NRPintronic = 0;
		this.NRPtranscriptomic = 0;
		this.NRPmergedTr = 0;
		this.NRPforHighestTranscript = 0;
	}
	
	
	public int getNRPwithinGenReg() {
		return NRPwithinGenReg;
	}
	public void addNRPwithinGenReg() {
		NRPwithinGenReg++;
	}
	
	public int getNRPwithinGenRegWitoutOtherGenReg() {
		return NRPwithinGenRegWitoutOtherGenReg;
	}
	public void addNRPwithinGenRegWitoutOtherGenReg() {
		NRPwithinGenRegWitoutOtherGenReg++;
	}
	
	public int getNRPintronic() {
		return NRPintronic;
	}
	public void addNRPintronic() {
		NRPintronic++;
	}
	
	public int getNRPtranscriptomic() {
		return NRPtranscriptomic;
	}
	public void addNRPtranscriptomic() {
		NRPtranscriptomic++;
	}
	
	public int getNRPmergedTr() {
		return NRPmergedTr;
	}
	public void addNRPmergedTr() {
		NRPmergedTr++;
	}
	
	public int getNRPforHighestTranscript() {
		return NRPforHighestTranscript;
	}
	public void addNRPforHighestTranscript() {
		NRPforHighestTranscript++;
	}
	
	
}
