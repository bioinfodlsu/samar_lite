use std::collections::{HashMap,HashSet};
use std::str;
use twox_hash::xxh3::hash64;

pub fn best_in_pair(p1: &Vec<Vec<(String,f64)>>, p2: &Vec<Vec<(String,f64)>>, frame1: usize, frame2: usize) -> String{
	let mut intersection: Vec<(String,f64)> = Vec::new();

	if frame1 < 10 && frame2 < 10{
		for x in &p1[frame1]{
			for y in &p2[frame2]{
				if x.0 == y.0{
					intersection.push(x.clone());
				}
			}
		}
	}
	else if frame2 < 10 {
		intersection = p2[frame2].clone();
	}
	else if frame1 < 10 {
		intersection = p1[frame1].clone();
	}
	let mut max = 0.0;
	let mut index = 99999;
	for (i,val) in intersection.iter().enumerate(){
		if val.1 > max {
			max = val.1;
			index = i;
		}
	}
	if index < 99999 {
		intersection[index].0.clone()
	}
	else{
		String::from("None")
	}
}

pub fn best_in_frame_max(frames: &Vec<Vec<(String,f64)>>) -> usize{
	let mut best = 10;
	let mut max = 0.0;

	for (i,info) in frames.iter().enumerate(){
		for j in info.iter(){
			if j.1 > max{
				max = j.1;
				best = i;
			}
		}
	}

	best
}

pub fn best_in_frame_sum(frames: &Vec<Vec<(String,f64)>>) -> usize{
	let mut best = 10;
	let mut max = 0.0;
	let mut sum = 0.0;

	for (i,info) in frames.iter().enumerate(){
		for j in info.iter(){
			sum+= j.1;
		}
		if sum > max{
			max = sum;
			best = i;
		}
		sum = 0.0;
	}

	best
}

pub fn mode(numbers: &[u64]) -> u64 {
    let mut occurrences = HashMap::new();

    for &value in numbers {
        *occurrences.entry(value).or_insert(0) += 1;
    }

    occurrences
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(val, _)| val)
        .expect("Cannot compute the mode of zero numbers")
}

pub fn b_search(s: &[u8], sa: &[usize], pat: &[u8]) -> (usize,usize) {
        
	//This section might be needed in the future if we need buckets 
	/*
	let sa = if pat.len() > 0 {
		&self.sa[self.get_bucket(pat)]
	} else {
		&self.sa[..]
	};
	*/
	b_search_mmp(s,sa,pat,0,sa.len())
}

pub fn b_search_mmp(s: &[u8], sa: &[usize], pat: &[u8], beg: usize, end: usize) -> (usize,usize) {
	let mut i = beg;
	let mut k = end;
	while i < k {
		let m = i + (k - i) / 2;
		if pat > &s[sa[m] as usize..] {
			i = m + 1;
		} else {
			k = m;
		}
	}

	let mut j = i;
	let mut k = end;
	while j < k {
		let m = j + (k - j) / 2;
		if s[sa[m] as usize..].starts_with(pat) {
			j = m + 1;
		} else {
			k = m;
		}
	}

	(i,j-1)
}

fn reduce_alph(seq:String) -> String{
	seq.chars()
        .map(|x| match x {
            'R' => 'K',
			'E' => 'K',
			'D' => 'K',
			'Q' => 'K',
			'N' => 'K',
			'L' => 'I',
			'V' => 'I',
			'T' => 'S',
			'A' => 'S',
            _ => x,
        })
        .collect()
}

// ToDo fix variable names, Match K to the reference index,argparse for threshold
pub fn palign(cat_str: &[u8], suff_arry: &[usize], pre_read:String, rs_maybe: &bio::data_structures::rank_select::RankSelect, ref_names: &HashMap<u64,String>, hash_table: &HashMap<u64,(u64,u64)>, bench: bool,kmer:usize, threshold: f64,  frame:u64) -> Vec<(String,f64)>{
	let mut trans: Vec<(String,f64)> = Vec::new();
	let mut cov_consensus: HashMap<String,f64> = HashMap::new(); 
	let mut x = 0;
	let almost_read = pre_read;
	let read = almost_read.as_bytes();
	let seqlen = read.len();
	// let align_time = Instant::now();
	// let kmer_time = Instant::now();
	//println!("frame {}: k {}: seq_len {}: READ {:?}",frame,kmer,seqlen,std::str::from_utf8(&read).unwrap());
	while x < (seqlen-kmer+1){
		let read_k = &read[x..(x+kmer)];
		//println!("Read: {:?}",read_k);
		// Go through suffix array and find a beg and end interval
		//let initial_bs = Instant::now();
		//let (beg,end) = b_search(cat_str, &suff_arry, read_k);

		//println!("{}",str_read_k);
		//println!("{:?}", hash_table.keys());
		let (beg,end) = if hash_table.contains_key(&hash64(&read_k)){
			let hash_k = &hash64(read_k);
			let mut sa_interval = hash_table.get(hash_k).unwrap().clone();
			let mut collision = 0;
			let mut hash_string = &cat_str[(suff_arry[sa_interval.0 as usize] as usize)..(suff_arry[sa_interval.0 as usize] as usize)+ kmer];
			
			//println!("INITIAL {}={}",str::from_utf8(read_k).unwrap(),str::from_utf8(hash_string).unwrap());
			while read_k != hash_string{
				//println!("IN LOOP {}={}",str::from_utf8(read_k).unwrap(),str::from_utf8(hash_string).unwrap());
			
				collision += 1;
				sa_interval = hash_table.get(&(hash_k + collision)).unwrap().clone();
				hash_string = &cat_str[(suff_arry[sa_interval.0 as usize] as usize)..(suff_arry[sa_interval.0 as usize] as usize)+ kmer];
			}
			// if collision != 0{
			// 	println!("{}",collision);
			// 	println!("IN LOOP {}={}",str::from_utf8(read_k).unwrap(),str::from_utf8(hash_string).unwrap());
			// }
			
			hash_table.get(&(hash_k + collision)).unwrap().clone()
		}
		else {
			(1,0)
		};
		// if bench{
		// 	println!("InitialBS: {:?}", initial_bs.elapsed());
		// }
		//println!("Align Initial Time: {:?} seconds", elapsed_align_initial);
		let mut mmp = 0;
		if beg <= end {
			//println!("Interval found: {}",read_k);
			//If an interval is found begin MMP
			let mut beg_p = beg as usize;
			let mut end_p = end as usize;
			
			//finding MMP
			//println!("SA index Match: Beg: {} & End: {}",beg,end);
			let mut done = false; 
			//let mmp_time = Instant::now();
			let mut ref_pro: HashSet<String> = HashSet::new();

			while !done{
				if (x+kmer+mmp+1) < read.len() {
					let temp_read = &read[x..(x+kmer+mmp)];
					let (temp_beg,temp_end) = b_search_mmp(cat_str, &suff_arry, temp_read,beg_p as usize,(end_p + 1) as usize);
					
					if temp_beg <= temp_end{
						beg_p = temp_beg;
						end_p = temp_end;
						mmp+=1;
					}
					else{
						done = true;
					}
				}
				else{
					done = true;
				}
			}
			// if bench{
			// 	println!("MMP: {:?}", mmp_time.elapsed());
			// }
			//println!("Align MMP Time: {:?} seconds", elapsed_align_mmp);
			//println!("MMP: {}",mmp);
			//println!("Beg prime: {} & End prime: {}", beg_p,end_p);
			//println!("{} {} {}",seqlen, kmer, mmp);
			//let add_consensus = Instant::now();
			//println!("mmp {}, x {}",mmp,x);
			//println!("PART OF READ {:?}", std::str::from_utf8(&read[x..(x+kmer+mmp)]).unwrap());
			for x in beg_p..(end_p+1){
				//consensus.push(HashSet::new().insert(rs_maybe.rank(suff_arry[x] as u64).unwrap()));
				//mode_consensus.push(rs_maybe.rank(suff_arry[x] as u64).unwrap());
				//if x >= beg_p as u64 && x <= end_p as u64{
				ref_pro.insert(ref_names[&rs_maybe.rank(suff_arry[x as usize]as u64).unwrap()].clone());
				/*
				}
				else
					*cov_consensus.entry(ref_names[&rs_maybe.rank(suff_arry[x as usize]as u64).unwrap()].clone()).or_insert(0.0) += (kmer) as f64 / seqlen as f64 * 100.0;
				}
				*/
			}
			//println!("Refs {:?}",ref_pro);
			// if bench{
			// 	println!("AddCon: {:?}",add_consensus.elapsed());
			// }
			for refs in ref_pro{
				*cov_consensus.entry(refs).or_insert(0.0) += (kmer + mmp) as f64;
			}
		}	
		if mmp > 0{
			x += kmer + mmp;
		}
		else{
			x += 1;
		}
	}
	// if bench{
	// 	println!("KmerPerRead: {:?}",kmer_time.elapsed());
	// }
	//Consesus machine 
	/*
	if !consensus.is_empty() {
		let (intersection, others) = consensus.split_at_mut(1);
		let intersection = &mut intersection[0];
		for other in others {
			intersection.retain(|e| other.contains(e));
		}
		//println!("{:?}",intersection);
		//Using Intersection
		
		if !intersection.is_empty(){
			*counts.entry(ref_names.get(intersection.iter().next().unwrap()).unwrap().to_string()).or_insert(0)  += 1;
		}
	*/
		// Using Coverage
		//else {
	if !cov_consensus.is_empty(){
		//println!{"Consensus: {:?}", cov_consensus};
		for (transcript, coverage) in cov_consensus.iter() {
			//println!("Transcript {}: Coverage {}",transcript,coverage);
			let cov = *coverage / seqlen as f64 * 100.0;
			//println!("Overall {}", cov);
			if cov > threshold {
				trans.push((transcript.clone(),cov));
			}
		}
		
	}
		//}
		/* Using Mode
		else {
			*counts.entry(ref_names.get(&mode(&mode_consensus)).unwrap().to_string()).or_insert(0)  += 1;
		}
		*/
	// if bench{
	// 	println!("AlignTime: {:?}", align_time.elapsed());
	// }
	trans
}
