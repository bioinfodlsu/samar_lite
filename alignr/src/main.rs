use std::path::Path;
use bio::io::fastq;
use std::str;
use bio::data_structures::rank_select::RankSelect;
use std::fs::File;
use bio::alphabets::dna::revcomp;
use protein_translate::translate;
use argparse::{ArgumentParser, Store};
use serde::{Deserialize};
use std::io::{BufReader,BufWriter,Write};
use std::collections::HashMap;

#[derive(Deserialize)]
struct Ref {
    cat_str: String,
    suff_arry: Vec<usize>,
    rank_select: RankSelect,
	ref_names: HashMap<u64,String>,
	hash_table: HashMap<u64, (u64,u64)>,
	k: usize
}

fn main() {
	let mut ref_file = String::from("ref_index_7.json");
	let mut query_file = String::from("sample_01_300000_interleaved.fastq");
	let mut aligns_file = String::from("output");

	let mut threshold:f64 = 30.0;

	{ // this block limits scope of borrows by ap.refer() method
		let mut ap = ArgumentParser::new();
		ap.set_description("Generate a Pseudo Alignment");
		ap.refer(&mut ref_file).add_argument("Reference Sequence", Store, "Reference Fasta file");
		ap.refer(&mut query_file).add_argument("Query Sequence", Store, "Query Interleaved Fastq file");
		ap.refer(&mut threshold).add_argument("Threshold", Store, "Coverage threshold (Must be Floating Point)");
		ap.refer(&mut aligns_file).add_argument("Alignment File", Store, "Output alignments of samar_lite");
		
		ap.parse_args_or_exit();
	}

	// let start = Instant::now();
	
	//For Reference
	// let start_ref = Instant::now();
	//TODO: PASS JSON LOCATION AND FILENAME FROM REFERENCe
	println!("BEGINNING IMPORT OF REF");
	let file = File::open(ref_file).unwrap();
    let reader = BufReader::new(file);

	let ref_struct:Ref = serde_json::from_reader(reader).unwrap();

	//For Query
	let path_q = Path::new(&query_file);
    let reader_q =  fastq::Reader::from_file(path_q);
	//let paths = fs::read_dir(r"D:/data/reads_20mil_1txpg_protCoding.rep3/query").unwrap();
	let mut records = reader_q.unwrap().records();
	let temp_cat_str = ref_struct.cat_str.as_bytes();
	let suff_arry = ref_struct.suff_arry;
	let rank_select = ref_struct.rank_select;
	let ref_names = ref_struct.ref_names;
	let hash_table = ref_struct.hash_table;
	let k = ref_struct.k;

	println!("DONE IMPORTING REF");

	let aligns = File::create(aligns_file).expect("unable to create output file");
    let mut b_aligns = BufWriter::new(aligns);

	println!("BEGINNING QUASI-MAPPING");
    while let (Some(Ok(pair1)),Some(Ok(pair2))) = (records.next(),records.next()){
		//println!("HERE IS READ: {} {}", pair1.id(),pair2.id());
		//println!("HERE IS READ: pair1{:?}\n pair2{:?}", pair1.seq(),pair2.seq());
		// let frame1_time = Instant::now();

		//Vector of 6 vectors, each element is a vector of tuples the gene name and the coverage 
		//println!("HERE IS READ1: {}",pair1.id());
		let p1_frame_con = vec![alignr::palign(temp_cat_str, &suff_arry, translate(&pair1.seq()), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&pair1.seq()[1..]), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&pair1.seq()[2..]), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())[1..]), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())[2..]), &rank_select, &ref_names, &hash_table,k,threshold)];
		
		// if bench{
		// 	println!("Frame1: {:?}",frame1_time.elapsed());
		// }
		//println!("FRAME1 CON \n{:?}",p1_frame_con);
		//let frame2_time = Instant::now();
		//println!("HERE IS READ2: {}",pair2.id());
		let p2_frame_con = vec![alignr::palign(temp_cat_str, &suff_arry, translate(&pair2.seq()), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&pair2.seq()[1..]), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&pair2.seq()[2..]), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq() ) ), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq() ) [1..]), &rank_select, &ref_names, &hash_table,k,threshold),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq() ) [2..]), &rank_select, &ref_names, &hash_table,k,threshold)];
		// if bench{
		// 	println!("Frame2: {:?}",frame2_time.elapsed());
		// }
		//IMPORTANT TODO: Utilize Coverage for Best in Frame (Max of Max and Max of Sum)
		//Notes: Pairing information can be used to improve choice

		//println!("FRAME2 CON \n{:?}",p2_frame_con);

		//let frame1_best_time = Instant::now();
		//let p1_frame = alignr::best_in_frame_sum(&p1_frame_con);
		// let p1_frame = alignr::best_in_frame_max(&p1_frame_con);
		// // if bench{
		// // 	println!("Frame1Best: {:?}",frame1_best_time.elapsed());
		// // }

		// //let frame2_best_time = Instant::now();
		// //let p2_frame = alignr::best_in_frame_sum(&p2_frame_con);
		// let p2_frame = alignr::best_in_frame_max(&p2_frame_con);
		// if bench{
		// 	println!("Frame2Best: {:?}",frame2_best_time.elapsed());
		// }
		
		let p1_frame = alignr::merge_frame(&p1_frame_con);
		let p2_frame = alignr::merge_frame(&p2_frame_con);
		

		// let best_pair = alignr::best_in_pair(&p1_frame_con, &p2_frame_con, p1_frame, p2_frame);
		
		// if best_pair != "None"{
		// 	*out.entry(best_pair).or_insert(0) += 1;
		// }

		//TODO: Segregate the pairs into p1 and p2 and into frames
		
		if !p1_frame.is_empty() {
			// let mut temp = Vec::new();
			// // if p1_frame.len() > 1{
			// // 	println!("MULTIFRAMEMAPPING!!!!");
			// // 	println!("{:?}",p1_frame);
			// // }

			// for x in p1_frame{
			// 	temp.push(p1_frame_con[x].clone());
			// }

			// // println!("OUTPUT: {:?}",temp);
			
			write!(b_aligns,"{}\t{:?}\n", pair1.id(), p1_frame).expect("unable to write");
		}

		if !p2_frame.is_empty() {
			// let mut temp = Vec::new();
			// // if p2_frame.len() > 1{
			// // 	println!("MULTIFRAMEMAPPING!!!!");
			// // 	println!("{:?}",p2_frame);
			// // }

			// for x in p2_frame{
			// 	temp.push(p2_frame_con[x].clone());
			// }

			// println!("OUTPUT: {:?}",temp);
			
			write!(b_aligns,"{}\t{:?}\n", pair2.id(), p2_frame).expect("unable to write");
		}
		
		//println!("{}",pair1.id());
		//println!("{}",pair2.id());
	} 

	// output_file.push_str(&k.to_string());
	// output_file.push_str(&threshold.to_string());
	// output_file.push_str(".count");
    // let f = File::create(output_file).expect("unable to create output file");
    // let mut f = BufWriter::new(f);

	// // Add Headers 22506112
	
	// write!(f,"Name\tNumReads\n").expect("unable to write");

	// for (key, value) in &out {
    //     write!(f,"{}\t{}\n", key, value).expect("unable to write");
    // }
}
