package bamfiles;

public class BAMStatistics {

	
	public void readBamFile(String bamFilePath){
		
		
		
	}
	
	
	
	public static void main(String[] args) {
		
		BAMFileReader br = new BAMFileReader("/home/proj/biosoft/praktikum/genprakt-ws16/assignment/a3/data/test_bams/star.bam");
		br.readBAMFile();
//		System.out.println(br.getCounter().toString());
		
	}
	
	
}
