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
	hash_table: HashMap<String, (u64,u64)>,
	k: usize
}

fn main() {
	let mut ref_file = String::from("test/data_suffarry_7.json");
	let mut query_file = String::from("data/sample_01_10000_interleaved.fastq");
	let mut output_file = String::from("output.counts");
	let mut aligns_file = String::from("output.aligns");

	let mut threshold:f64 = 30.0;
	let mut bench = false;

	{ // this block limits scope of borrows by ap.refer() method
		let mut ap = ArgumentParser::new();
		ap.set_description("Generate a Pseudo Alignment");
		ap.refer(&mut ref_file).add_argument("Reference Sequence", Store, "Reference Fasta file");
		ap.refer(&mut query_file).add_argument("Query Sequence", Store, "Query Interleaved Fastq file");
		ap.refer(&mut threshold).add_argument("Threshold", Store, "Coverage threshold (Must be Floating Point)");
		ap.refer(&mut aligns_file).add_argument("Alignment File", Store, "Output alignments of samar_lite");
		ap.refer(&mut output_file).add_argument("Count File", Store, "Output counts of samar_lite");
		
		ap.parse_args_or_exit();
	}

	// let start = Instant::now();
	
	//For Reference
	// let start_ref = Instant::now();
	//TODO: PASS JSON LOCATION AND FILENAME FROM REFERENCe
	let file = File::open(ref_file).unwrap();
    let reader = BufReader::new(file);

	let ref_struct:Ref = serde_json::from_reader(reader).unwrap();

	// let elapsed_ref = start_ref.elapsed();
	/* Print Checks
	println!("Number of reads: {}", nb_reads);
	println!("Number of bases: {}", nb_bases);
	println!("Concatenation of all Reads: \n{}", cat_str);
	println!("Len of SA: {}",suff_arry.len());
	println!("Index of SA: {}",suff_arry[126]);
	println!("Index of SA: {}",suff_arry[127]);
	println!("Rank of element 1300 {:?}",rs_maybe.rank(suff_arry[126]as u64));
	println!("Rank of element 1300 {:?}",rs_maybe.rank(suff_arry[127]as u64));
	*/

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

	let mut out: HashMap<String,u64> = HashMap::new();	
	aligns_file.push_str(&k.to_string());
	aligns_file.push_str(&threshold.to_string());
	let aligns = File::create(aligns_file).expect("unable to create output file");
    let mut b_aligns = BufWriter::new(aligns);

    while let (Some(Ok(pair1)),Some(Ok(pair2))) = (records.next(),records.next()){
		//println!("HERE IS READ: {} {}", pair1.id(),pair2.id());
		//println!("HERE IS READ: pair1{:?}\n pair2{:?}", pair1.seq(),pair2.seq());
		// let frame1_time = Instant::now();

		//Vector of 6 vectors, each element is a vector of tuples the gene name and the coverage 
		//println!("HERE IS READ1: {}",pair1.id());
		let p1_frame_con = vec![alignr::palign(temp_cat_str, &suff_arry, translate(&pair1.seq()), &rank_select, &ref_names, &hash_table,bench,k,threshold,1),
								alignr::palign(temp_cat_str, &suff_arry, translate(&pair1.seq()[1..]), &rank_select, &ref_names, &hash_table,bench,k,threshold,2),
								alignr::palign(temp_cat_str, &suff_arry, translate(&pair1.seq()[2..]), &rank_select, &ref_names, &hash_table,bench,k,threshold,3),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())), &rank_select, &ref_names, &hash_table,bench,k,threshold,4),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())[1..]), &rank_select, &ref_names, &hash_table,bench,k,threshold,5),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair1.seq())[2..]), &rank_select, &ref_names, &hash_table,bench,k,threshold,6)];
		
		// if bench{
		// 	println!("Frame1: {:?}",frame1_time.elapsed());
		// }
		//println!("FRAME1 CON \n{:?}",p1_frame_con);
		//let frame2_time = Instant::now();
		//println!("HERE IS READ2: {}",pair2.id());
		let p2_frame_con = vec![alignr::palign(temp_cat_str, &suff_arry, translate(&pair2.seq()), &rank_select, &ref_names, &hash_table,bench,k,threshold,1),
								alignr::palign(temp_cat_str, &suff_arry, translate(&pair2.seq()[1..]), &rank_select, &ref_names, &hash_table,bench,k,threshold,2),
								alignr::palign(temp_cat_str, &suff_arry, translate(&pair2.seq()[2..]), &rank_select, &ref_names, &hash_table,bench,k,threshold,3),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq() ) ), &rank_select, &ref_names, &hash_table,bench,k,threshold,4),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq() ) [1..]), &rank_select, &ref_names, &hash_table,bench,k,threshold,5),
								alignr::palign(temp_cat_str, &suff_arry, translate(&revcomp(pair2.seq() ) [2..]), &rank_select, &ref_names, &hash_table,bench,k,threshold,6)];
		// if bench{
		// 	println!("Frame2: {:?}",frame2_time.elapsed());
		// }
		//IMPORTANT TODO: Utilize Coverage for Best in Frame (Max of Max and Max of Sum)
		//Notes: Pairing information can be used to improve choice

		//println!("FRAME2 CON \n{:?}",p2_frame_con);

		//let frame1_best_time = Instant::now();
		//let p1_frame = alignr::best_in_frame_sum(&p1_frame_con);
		let p1_frame = alignr::best_in_frame_max(&p1_frame_con);
		// if bench{
		// 	println!("Frame1Best: {:?}",frame1_best_time.elapsed());
		// }

		//let frame2_best_time = Instant::now();
		//let p2_frame = alignr::best_in_frame_sum(&p2_frame_con);
		let p2_frame = alignr::best_in_frame_max(&p2_frame_con);
		// if bench{
		// 	println!("Frame2Best: {:?}",frame2_best_time.elapsed());
		// }

		let best_pair = alignr::best_in_pair(&p1_frame_con, &p2_frame_con, p1_frame, p2_frame);
		
		if best_pair != "None"{
			*out.entry(best_pair).or_insert(0) += 1;
		}
		/*
		//println!("{:?}",p1_frame);
		//println!("{:?}",p2_frame);
		let best_pair_time = Instant::now();
		;
		if bench{
		println!("BestPair: {:?}",best_pair_time.elapsed());
		}
		println!("{} Aligns to {:?}",pair1.id(), alignr::best_in_pair(&p1_frame_con, &p2_frame_con, p1_frame, p2_frame));
		*/

		// let temp_consesus = Instant::now();
		// let mut last1:Vec<(String,f64)> = Vec::new();
		// let mut last2:Vec<(String,f64)> = Vec::new();

		// for mut i in p1_frame_con[p1_frame]{
		// 	last1.append(&mut i);
		// }

		// for mut i in p2_frame_con[p2_frame]{
		// 	if !i.is_empty(){
		// 		last2.append(&mut i);
		// 	}
		// }

		// if bench{
		// 	println!("TempCon: {:?}",temp_consesus.elapsed());
		// }

		//TODO: Segregate the pairs into p1 and p2 and into frames
		
		if p1_frame < 10 {
			write!(b_aligns,"{}\t{:?}\n", pair1.id(), p1_frame_con[p1_frame]).expect("unable to write");
		}

		if p2_frame < 10 {
			write!(b_aligns,"{}\t{:?}\n", pair2.id(), p2_frame_con[p2_frame]).expect("unable to write");
		}
		
		//println!("{}",pair1.id());
		//println!("{}",pair2.id());
	} 

	//println!("Count: {}", count);

	/*
	let elapsed_query = start_query.elapsed();
	//println!("{:?}", read_alignments);
	let elapsed = start.elapsed();
	if bench{
		println!("Total: {:?}", elapsed);
		println!("Query: {:?}", elapsed_query);
		println!("Ref: {:?}", elapsed_ref);
	}
	*/
	output_file.push_str(&k.to_string());
	output_file.push_str(&threshold.to_string());
    let f = File::create(output_file).expect("unable to create output file");
    let mut f = BufWriter::new(f);

	// Add Headers 22506112
	
	write!(f,"Name\tNumReads\n").expect("unable to write");

	for (key, value) in &out {
        write!(f,"{}\t{}\n", key, value).expect("unable to write");
    }
}