use std::str;
use std::fs::File;
use serde::{Serialize};
use std::path::Path;
use bio::io::fasta;
use bio::data_structures::qgram_index;
use std::collections::HashMap;
use bv::BitVec;
use bio::data_structures::rank_select::RankSelect;
use bio::alphabets::Alphabet;
use bio::alphabets;
use bio::data_structures::suffix_array::suffix_array;
use argparse::{ArgumentParser, StoreTrue, Store};

// Reference structure
#[derive(Serialize)]
struct Ref {
    cat_str: String,
    //qgram: qgram_index::QGramIndex,
	suff_arry: Vec<usize>,
    rank_select: RankSelect,
	ref_names: HashMap<u64,String>,
	hash_table: HashMap<String,(u64,u64)>,
	k: usize
}

fn main() {
	/*To Do: 
		Add ArgParse to Input the reference
	*/

	// Default values for reference file path and k-mer size
	let mut ref_path = String::from("ref.fasta");
	let mut k:usize = 5;

	{ // this block limits scope of borrows by ap.refer() method
		let mut ap = ArgumentParser::new();
		ap.set_description("Reference Creator for Samar-lite");
		ap.refer(&mut ref_path).add_argument("Reference", Store,"Reference Fasta File");
		ap.refer(&mut k).add_argument("kmer", Store,"Kmer for alignment");
		ap.parse_args_or_exit();
	}
	/*
		COMMENT FOR ^^:
		It doesn't have an arguement parser for the
		path to output of the json file containing
		the suffix array. Instead the name of the file
		and its path (./) are directly passed as a string
		literal on line 122 (serde_json::to_writer() function call)
	*/

	let path_r = Path::new(&ref_path);
    let reader = fasta::Reader::from_file(path_r).expect("File not found");

	let mut ref_names: HashMap<u64,String> = HashMap::new();
	let mut cat_str = String::new();
	let mut rank_bitvec: BitVec<u8> = BitVec::new();
	
	let q = 5 as u32; // unused variables q and alph
	let alph = Alphabet::new(&b"ABCDEFGHIJKLMNOPQRSTUVWXYZ$*abcdefghijklmnopqrstuvwxyz"[..]);

	// Creates Concatenated string and Rank Bit Vector and Rank Names Map
	let mut nb_reads = 0;

	for result in reader.records() {
		let record = result.expect("Error during fasta record parsing");
		
		cat_str.push_str(str::from_utf8(record.seq()).expect("Error during fasta record parsing"));
		cat_str.push_str("$");
		
		for _x in 0..record.seq().len(){
			rank_bitvec.push(false);
		}
		
		rank_bitvec.push(true);
		ref_names.insert(nb_reads,record.id().to_string());
		nb_reads += 1;
	}

	let suff_arry = suffix_array(&cat_str.as_bytes());
	
	// K-mer to Suffix Array Interval Hash Table Construction

	let n = cat_str.len();
	let mut hash_table: HashMap<String,(u64,u64)> = HashMap::new();
	let mut prev_str = "";
	
	let mut j = "";
	let mut right: u64 = 0;
	let mut left: u64 = 0;

	for x in 0..n{
		if suff_arry[x] <= (n - k){
			//println!("{} {}", suff_arry[x],cat_str.len());
			//println!("{}",cat_str[suff_arry[x]..suff_arry[x]+k-1].to_owned());

			if &cat_str[suff_arry[x]..suff_arry[x]+k] != prev_str{
				if j != ""{
					right = (x - 1) as u64;
					hash_table.insert(j.to_owned(),(left,right));
				}
				j = &cat_str[suff_arry[x]..suff_arry[x]+k];
				prev_str = &j;
				left = x as u64;
			}
		}
	}
	hash_table.insert(j.to_owned(),(left,(n-1) as u64));
	
	let rank_select = RankSelect::new(rank_bitvec,1); 
	
	let test = Ref{
		cat_str,
		//qgram,
		suff_arry,
		rank_select,
		ref_names,
		hash_table,
		k
	};
	
	serde_json::to_writer(&File::create("data_suffarry.json").unwrap(), &test);
	// ^ Maybe wrap in a try-catch block for safety since it returns a Result type.
	
	println!("Done");
}
