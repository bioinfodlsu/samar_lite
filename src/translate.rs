use std::path::Path;
use std::fs;
use bio::io::fasta;
use bio::io::fastq;
use bio::alphabets::rna::revcomp;
use protein_translate::translate;
use std::str;

pub fn p_trans(path_r :&Path) {
	let file_name = path_r.file_name().unwrap().to_str().unwrap();
	let name_only = &file_name[..file_name.len()-6];
	
    //let path_r = Path::new("data/query/reads_1_1.fastq");
	/*
	let folder_name = "translated/".to_owned() + name_only;
	fs::create_dir(folder_name.clone());
	
	let path1 = "".to_owned() + &folder_name + "/" + name_only +  "_frame1.fasta";
	let path2 = "".to_owned() + &folder_name + "/" + name_only +  "_frame2.fasta";
	let path3 = "".to_owned() + &folder_name + "/" + name_only +  "_frame3.fasta";
	let path4 = "".to_owned() + &folder_name + "/" + name_only +  "_frame4.fasta";
	let path5 = "".to_owned() + &folder_name + "/" + name_only +  "_frame5.fasta";
	let path6 = "".to_owned() + &folder_name + "/" + name_only +  "_frame6.fasta";
	*/
	let path1 = "translated/".to_owned()  + name_only +  "_frame1.fasta";
	let path2 = "translated/".to_owned()  + name_only +  "_frame2.fasta";
	let path3 = "translated/".to_owned()  + name_only +  "_frame3.fasta";
	let path4 = "translated/".to_owned()  + name_only +  "_frame4.fasta";
	let path5 = "translated/".to_owned()  + name_only +  "_frame5.fasta";
	let path6 = "translated/".to_owned()  + name_only +  "_frame6.fasta";
	
	let path_w_f1 = Path::new(&path1);
	let path_w_f2 = Path::new(&path2);
	let path_w_f3 = Path::new(&path3);
	let path_w_f4 = Path::new(&path4);
	let path_w_f5 = Path::new(&path5);
	let path_w_f6 = Path::new(&path6);
	
	let file_f1 = fs::File::create(path_w_f1).unwrap();
	let file_f2 = fs::File::create(path_w_f2).unwrap();
	let file_f3 = fs::File::create(path_w_f3).unwrap();
	let file_f4 = fs::File::create(path_w_f4).unwrap();
	let file_f5 = fs::File::create(path_w_f5).unwrap();
	let file_f6 = fs::File::create(path_w_f6).unwrap();
	
    let reader = fastq::Reader::from_file(path_r).expect("File not found");
	let mut writer_f1 = fasta::Writer::new(file_f1);
	let mut writer_f2 = fasta::Writer::new(file_f2);
	let mut writer_f3 = fasta::Writer::new(file_f3);
	let mut writer_f4 = fasta::Writer::new(file_f4);
	let mut writer_f5 = fasta::Writer::new(file_f5);
	let mut writer_f6 = fasta::Writer::new(file_f6);
	
	let mut nb_reads = 0;
	let mut nb_bases = 0;
	
	for result in reader.records() {
		let record = result.expect("Error during fastq record parsing");
		
		let frame_1 = "frame_1_".to_owned() + record.id();
		let frame_2 = "frame_2_".to_owned() + record.id();
		let frame_3 = "frame_3_".to_owned() + record.id();
		let frame_4 = "frame_4_".to_owned() + record.id();
		let frame_5 = "frame_5_".to_owned() + record.id();
		let frame_6 = "frame_6_".to_owned() + record.id();
		
		writer_f1.write(&frame_1, None, translate(record.seq()).as_bytes()).ok().expect("Error writing record.");
		writer_f2.write(&frame_2, None, translate(&record.seq()[1..]).as_bytes()).ok().expect("Error writing record.");
		writer_f3.write(&frame_3, None, translate(&record.seq()[2..]).as_bytes()).ok().expect("Error writing record.");
		writer_f4.write(&frame_4, None, translate(&revcomp(record.seq())).as_bytes()).ok().expect("Error writing record.");
		writer_f5.write(&frame_5, None, translate(&revcomp(record.seq())[1..]).as_bytes()).ok().expect("Error writing record.");
		writer_f6.write(&frame_6, None, translate(&revcomp(record.seq())[2..]).as_bytes()).ok().expect("Error writing record.");
		
		nb_reads += 1;
		nb_bases += record.seq().len();
	}
	
	//println!("Number of reads: {}", nb_reads);
	//println!("Number of bases: {}", nb_bases);
}