use std::str;
use std::fs::File;
use serde::{Serialize};
use std::path::Path;
use bio::io::fasta;
use std::collections::HashMap;
use bv::BitVec;
use bio::data_structures::rank_select::RankSelect;
use bio::alphabets::Alphabet;
use bio::data_structures::suffix_array::suffix_array;
use argparse::{ArgumentParser, Store, StoreTrue};
use std::hash::BuildHasherDefault;
use twox_hash::xxh3::hash64;

#[derive(Serialize)]
struct Ref {
    cat_str: String,
    //qgram: qgram_index::QGramIndex,
	suff_arry: Vec<usize>,
    rank_select: RankSelect,
	ref_names: HashMap<u64,String>,
	hash_table: HashMap<u64, (u64,u64)>,
	k: usize
	// reduced: bool
}

// fn reduce_alph(seq:String) -> String{
// 	seq.chars().map(|x| match x {
// 		'R' => 'K',
// 		'E' => 'K',
// 		'D' => 'K',
// 		'Q' => 'K',
// 		'N' => 'K',
// 		'L' => 'I',
// 		'V' => 'I',
// 		'T' => 'S',
// 		'A' => 'S',
// 		_ => x,
// 	})
// 	.collect()
// }


fn main() {
	/*To Do: 
		Add ArgParse to Input the reference
	*/
	let mut ref_path = String::from("ref.fasta");
	let mut k:usize = 5;
	// let mut reduced = false;
	let mut output = String::from("ref_index_");
	output.push_str(&k.to_string());
	output.push_str(".json");

	{ // this block limits scope of borrows by ap.refer() method
	let mut ap = ArgumentParser::new();
	ap.set_description("Reference Creator for Samar-lite");
	ap.refer(&mut ref_path).add_argument("Reference", Store,"Reference Fasta File");
	ap.refer(&mut k).add_argument("kmer", Store,"Kmer for alignment");
	ap.refer(&mut output).add_argument("output", Store, "Specify output file");
	ap.parse_args_or_exit();
	}

	let path_r = Path::new(&ref_path);
    let reader = fasta::Reader::from_file(path_r).expect("File not found");
	let mut ref_names: HashMap<u64,String> = HashMap::new();
	let mut cat_str = String::new();
	let mut rank_bitvec: BitVec<u8> = BitVec::new();
	
	// if reduced{
	// 	println!("Using Reduced Alph");	
	// }

	//Creates Concatenated string and Rank Bit Vector and Rank Names Map
	let mut nb_reads = 0;

	for result in reader.records() {
		let record = result.expect("Error during fasta record parsing");
		
		// if reduced{
		// 	cat_str.push_str(&reduce_alph(String::from_utf8_lossy(record.seq()).into_owned()));
		// }
		// else {
		cat_str.push_str(str::from_utf8(record.seq()).expect("Error during fasta record parsing"));
		// }
		cat_str.push_str("$");
		
		for _x in 0..record.seq().len(){
			rank_bitvec.push(false);
		}
		
		rank_bitvec.push(true);
		ref_names.insert(nb_reads,record.id().to_string());
		nb_reads += 1;
	}

	println!("DoNE WITH RANK");

	let suff_arry = suffix_array(&cat_str.as_bytes());
	
	println!("DONE WITH SA");

	//Kmer to Suffix Array Interval Hash Table Construction

	let n = cat_str.len();
	let mut hash_table: HashMap<u64, (u64,u64)> = HashMap::new();
	let mut prev_str = "";
	
	let mut j = u64::MAX;
	let mut right: u64 = 0;
	let mut left: u64 = 0;

	for x in 0..n{
		if suff_arry[x] <= (n - k){
			//println!("{} {}", suff_arry[x],cat_str.len());
			//println!("{}",cat_str[suff_arry[x]..suff_arry[x]+k-1].to_owned());

			if &cat_str[suff_arry[x]..suff_arry[x]+k] != prev_str{
				if j != u64::MAX{
					right = (x - 1) as u64;
					hash_table.insert(j,(left,right));
				}
				j = hash64(&cat_str[suff_arry[x]..suff_arry[x]+k].as_bytes());
				prev_str = &cat_str[suff_arry[x]..suff_arry[x]+k];

				while hash_table.contains_key(&j){
					j+=1; 
				}
				
				left = x as u64;
			}
		}
	}
	hash_table.insert(j,(left,(n-1) as u64));

	let rank_select = RankSelect::new(rank_bitvec,1); 
	
	let test = Ref{
		cat_str,
		//qgram,
		suff_arry,
		rank_select,
		ref_names,
		hash_table,
		k
		// reduced
	};
		
	serde_json::to_writer(&File::create(output).unwrap(), &test);
	
	println!("Done");
}
