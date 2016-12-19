package assignment_3;

import bamfiles.BAMFileReader;
import bamfiles.Counter;

public class Task_1 {

	public static void main(String[] args) {
		new BAMFileReader("G:/contextmap.bam", new Counter(), null).readBAMFile();
	}
}