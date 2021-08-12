use std::collections::HashMap;
use std::time::Instant;
use std::str;

pub fn best_in_pair(p1: &Vec<Vec<(String,f64)>>, p2: &Vec<Vec<(String,f64)>>, frame1: usize, frame2: usize) -> Vec<String>{
	let mut matches: Vec<String> = Vec::new();

	if frame1 < 10 && frame2 < 10{
		for x in &p1[frame1]{
			for y in &p2[frame2]{
				if x.0 == y.0{
					matches.push(x.0.clone());
				}
			}
		}
	}
	else if frame2 < 10 {
		for x in &p2[frame2]{
			matches.push(x.0.clone());
		}
	}
	else if frame1 < 10 {
		for x in &p1[frame1]{
			matches.push(x.0.clone());
		}
	}

	matches
}

pub fn best_in_frame(frames: &Vec<Vec<(String,f64)>>) -> usize{
	let mut best = 10;
	let mut max = 0;

	for (i,info) in frames.iter().enumerate(){
		if info.len() > max{
			max = info.len();
			best = i;
		}
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



pub fn palign(cat_str: &[u8], suff_arry: &[usize], read:&[u8], rs_maybe: &bio::data_structures::rank_select::RankSelect, ref_names: &HashMap<u64,String>, hash_table: &HashMap<String,(u64,u64)>, bench: bool) -> Vec<(String,f64)>{
	let mut trans: Vec<(String,f64)> = Vec::new();
	let kmer = 5;
	let mut cov_consensus: HashMap<String,f64> = HashMap::new(); 
	let threshold:f64 = 1.0;
	let mut x = 0;
	let seqlen = read.len();
	let align_time = Instant::now();
	let kmer_time = Instant::now();
	while x < (seqlen-kmer+1){
		let read_k = &read[x..(x+kmer)];
		//println!("{:?}",read_k);
		// Go through suffix array and find a beg and end interval
		let initial_bs = Instant::now();
		//let (beg,end) = b_search(cat_str, &suff_arry, read_k);
		let str_read_k = str::from_utf8(read_k).unwrap();

		//println!("{}",str_read_k);
		//println!("{:?}", hash_table.keys());
		let (beg,end) = if hash_table.contains_key(str_read_k){
			hash_table.get(str_read_k).unwrap().clone()
		}
		else {
			(1,0)
		};

		//println!("{:?} {:?}", beg, end);
		if bench{
			println!("InitialBS: {:?}", initial_bs.elapsed());
		}
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
			let mmp_time = Instant::now();
			while !done{
				if (x+kmer+mmp+1) < read.len() {
					let temp_read = &read[x..(x+kmer+mmp+1)];
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
			if bench{
				println!("MMP: {:?}", mmp_time.elapsed());
			}
			//println!("Align MMP Time: {:?} seconds", elapsed_align_mmp);
			//println!("MMP: {}",mmp);
			//println!("Beg prime: {} & End prime: {}", beg_p,end_p);
			//println!("{} {} {}",seqlen, kmer, mmp);
			let add_consensus = Instant::now();
			for x in beg_p..(end_p+1){
				//consensus.push(HashSet::new().insert(rs_maybe.rank(suff_arry[x] as u64).unwrap()));
				//mode_consensus.push(rs_maybe.rank(suff_arry[x] as u64).unwrap());
				*cov_consensus.entry(ref_names[&rs_maybe.rank(suff_arry[x] as u64).unwrap()].clone()).or_insert((kmer + mmp) as f64 / seqlen as f64 * 100.0) += 0.0 ;
			}
			if bench{
				println!("AddCon: {:?}",add_consensus.elapsed());
			}
		}	
		x += 1+mmp;
	}
	if bench{
		println!("KmerPerRead: {:?}",kmer_time.elapsed());
	}
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
		let consensus_time = Instant::now();
		for (transcript, coverage) in cov_consensus.iter() {
			if *coverage > threshold {
				trans.push((transcript.clone(),*coverage));
			}
		}
		//
		if bench{
			println!("Consensus: {:?}",consensus_time.elapsed());
		}
	}
		//}
		/* Using Mode
		else {
			*counts.entry(ref_names.get(&mode(&mode_consensus)).unwrap().to_string()).or_insert(0)  += 1;
		}
		*/
	if bench{
		println!("AlignTime: {:?}", align_time.elapsed());
	}
	trans
}
