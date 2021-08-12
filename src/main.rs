use std::time::Instant;
use std::path::Path;
use bio::io::fastq;
use std::str;
use bio::data_structures::rank_select::RankSelect;
use std::fs::File;
use bio::alphabets::rna::revcomp;
use protein_translate::translate;
use argparse::{ArgumentParser, StoreTrue, Store};
use serde::{Deserialize};
use std::io::BufReader;
use std::collections::HashMap;

#[derive(Deserialize)]
struct Ref {
    cat_str: String,
    suff_arry: Vec<usize>,
    rs_maybe: RankSelect,
	ref_names: HashMap<u64,String>,
	hash_table: HashMap<String, (u64,u64)>
}

fn main() {

	let mut ref_file = String::from("data/ref.fasta");
	let mut query_file = String::from("data/sample_01_300000_interleaved.fastq");
	// { // this block limits scope of borrows by ap.refer() method
	// let mut ap = ArgumentParser::new();
	// ap.set_description("Greet somebody.");
	// ap.refer(&mut ref_file)
	// .add_argument("Reference Sequence", Store, "Reference Fasta file");
	// ap.refer(&mut query_file)
	// .add_argument("Query Sequence", Store, "Query Interleaved Fastq file");
	// ap.parse_args_or_exit();
	// }
	let bench = false;

	let start = Instant::now();
	
	//For Reference
	let start_ref = Instant::now();
	
	let file = File::open("data/data_suffarry.json").unwrap();
    let reader = BufReader::new(file);

	let ref_struct:Ref = serde_json::from_reader(reader).unwrap();

	let elapsed_ref = start_ref.elapsed();
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
	let start_query = Instant::now();
	let path_q = Path::new(&query_file);
    let reader_q =  fastq::Reader::from_file(path_q);
	//let paths = fs::read_dir(r"D:/data/reads_20mil_1txpg_protCoding.rep3/query").unwrap();
	let mut records = reader_q.unwrap().records();
	
    while let (Some(Ok(pair1)),Some(Ok(pair2))) = (records.next(),records.next()){
		//println!("HERE IS READ: {} {}", pair1.id(),pair2.id());

		let temp_cat_str = ref_struct.cat_str.as_bytes();

		let frame1_time = Instant::now();

		let p1_frame_con = vec![alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(pair1.seq()).as_bytes(), &ref_struct.rs_maybe, &ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&pair1.seq()[1..]).as_bytes(), &ref_struct.rs_maybe, &ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&pair1.seq()[2..]).as_bytes(), &ref_struct.rs_maybe, &ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&revcomp(pair1.seq())).as_bytes(), &ref_struct.rs_maybe, &ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&revcomp(pair1.seq())[1..]).as_bytes(), &ref_struct.rs_maybe, &ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&revcomp(pair1.seq())[2..]).as_bytes(), &ref_struct.rs_maybe, &ref_struct.ref_names, &ref_struct.hash_table,bench)];
		
		if bench{
			println!("Frame1: {:?}",frame1_time.elapsed());
		}
		//println!("{:?}",p1_frame_con);
		let frame2_time = Instant::now();

		let p2_frame_con = vec![alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(pair2.seq()).as_bytes(), &ref_struct.rs_maybe,&ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&pair2.seq()[1..]).as_bytes(), &ref_struct.rs_maybe,&ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&pair2.seq()[2..]).as_bytes(), &ref_struct.rs_maybe,&ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&revcomp(pair2.seq())).as_bytes(), &ref_struct.rs_maybe,&ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&revcomp(pair2.seq())[1..]).as_bytes(), &ref_struct.rs_maybe,&ref_struct.ref_names, &ref_struct.hash_table,bench),
								alignr::palign(temp_cat_str, &ref_struct.suff_arry, translate(&revcomp(pair2.seq())[2..]).as_bytes(), &ref_struct.rs_maybe,&ref_struct.ref_names, &ref_struct.hash_table,bench)];
		if bench{
			println!("Frame2: {:?}",frame2_time.elapsed());
		}

		/*
		println!("{:?}",p2_frame_con);

		let frame1_best_time = Instant::now();
		let p1_frame = alignr::best_in_frame(&p1_frame_con);
		if bench{
			println!("Frame1Best: {:?}",frame1_best_time.elapsed());
		}

		let frame2_best_time = Instant::now();
		let p2_frame = alignr::best_in_frame(&p2_frame_con);
		if bench{
			println!("Frame2Best: {:?}",frame2_best_time.elapsed());
		}

		//println!("{:?}",p1_frame);
		//println!("{:?}",p2_frame);
		let best_pair_time = Instant::now();
		alignr::best_in_pair(&p1_frame_con, &p2_frame_con, p1_frame, p2_frame);
		if bench{
		println!("BestPair: {:?}",best_pair_time.elapsed());
		}
		println!("{} Aligns to {:?}",pair1.id(), alignr::best_in_pair(&p1_frame_con, &p2_frame_con, p1_frame, p2_frame));
		*/

		let temp_consesus = Instant::now();
		let mut last = Vec::new();

		for mut i in p1_frame_con{
			if !i.is_empty(){
				last.append(&mut i);
			}
		}

		for mut i in p2_frame_con{
			if !i.is_empty(){
				last.append(&mut i);
			}
		}

		if bench{
			println!("TempCon: {:?}",temp_consesus.elapsed());
		}

		//println!("{} Aligns to {:?}",pair1.id(), last);
	} 
	let elapsed_query = start_query.elapsed();
	//println!("{:?}", read_alignments);
	let elapsed = start.elapsed();
	if bench{
		println!("Total: {:?}", elapsed);
		println!("Query: {:?}", elapsed_query);
		println!("Ref: {:?}", elapsed_ref);
	}
}