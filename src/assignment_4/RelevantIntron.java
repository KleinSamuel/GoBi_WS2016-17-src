package assignment_4;

import java.util.ArrayList;

import util.Interval;

public class RelevantIntron extends Interval{

	private String geneID;
	private String chromosomeID;
	private int start;
	private int stop;
	private boolean isOnNegativeStrand;
	private ArrayList<RelevantIntron> overlappingIntrons;
	boolean toRemove = true;
	double[] values;
	int pointerInValues;
	
	public RelevantIntron(String geneID, String chromosomeId, int start, int stop, boolean isOnNegativeStrand){
		super(start, stop);
		this.setGeneID(geneID);
		this.setChromosomeID(chromosomeId);
		this.setStart(start);
		this.setStop(stop);
		this.setOnNegativeStrand(isOnNegativeStrand);
		this.setOverlappingIntrons(new ArrayList<>());
	}
	
	public boolean isTheSame(RelevantIntron i){
		return (isOnSameGene(i.getGeneID()) && sameRegion(i.getStart(), i.getStop()) && isOnSameStrand(i.isOnNegativeStrand()));
	}
	
	public boolean isOnSameGene(String geneId){
		return this.geneID.equals(geneId);
	}
	
	public boolean sameRegion(int start, int stop){
		return (this.start == start && this.stop == stop);
	}
	
	public boolean isOnSameStrand(boolean isNegative){
		return (this.isOnNegativeStrand == isNegative);
	}
	
	public boolean isOnSameChromosome(String chrID){
		return this.chromosomeID.equals(chrID);
	}
	
	public boolean overlaps(RelevantIntron i){
		if(i.getStart() < this.start && i.getStop() >= this.start){
			return true;
		}else if(i.getStart() < this.stop && i.getStop() > this.getStop()){
			return true;
		}else if(i.getStart() < this.start && i.getStop() > this.stop){
			return true;
		}else if(this.start < i.getStart() && this.stop > i.getStop()){
			return true;
		}
		return false;
	}

	public String getGeneID() {
		return geneID;
	}

	public void setGeneID(String geneID) {
		this.geneID = geneID;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getStop() {
		return stop;
	}

	public void setStop(int stop) {
		this.stop = stop;
	}

	public boolean isOnNegativeStrand() {
		return isOnNegativeStrand;
	}

	public void setOnNegativeStrand(boolean isOnNegativeStrand) {
		this.isOnNegativeStrand = isOnNegativeStrand;
	}

	public ArrayList<RelevantIntron> getOverlappingIntrons() {
		return overlappingIntrons;
	}

	public void setOverlappingIntrons(ArrayList<RelevantIntron> overlappingIntrons) {
		this.overlappingIntrons = overlappingIntrons;
	}
	
	public void addOverlappingIntron(RelevantIntron i){
		this.overlappingIntrons.add(i);
	}

	public String getChromosomeID() {
		return chromosomeID;
	}

	public void setChromosomeID(String chromosomeID) {
		this.chromosomeID = chromosomeID;
	}
	
}
