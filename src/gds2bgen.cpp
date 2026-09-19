// ===========================================================
//
// gds2bgen.cpp: format conversion between GDS and BGEN
//
// Copyright (C) 2018-2025    Xiuwen Zheng (zhengxwen@gmail.com)
//
// gds2bgen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// gds2bgen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with gds2bgen.
// If not, see <http://www.gnu.org/licenses/>.


#include <R_GDS_CPP.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define class xclass
#define private xprivate
#include <R_ext/Connections.h>
#undef class
#undef private

#include "genfile/bgen/bgen.hpp"
#include <fstream>


using namespace std;
using namespace CoreArray;
using namespace genfile;


// defined in bgen_v1.1.8/wscript
static const char *bgen_version = "bgen_lib_v1.1.8";

// The ProbSetter structure is from bgen_to_vcf.cpp
// ProbSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comment in bgen.hpp above
// bgen::read_genotype_data_block(), or the bgen wiki for a description of the API.
// The purpose of this object is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
struct ProbSetter {
	typedef std::vector< std::vector< double > > Data ;
	ProbSetter( Data* result ):
		m_result( result ),
		m_sample_i(0)
	{}
		
	// Called once allowing us to set storage.
	void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
		m_result->clear() ;
		m_result->resize( number_of_samples ) ;
	}
	
	// If present with this signature, called once after initialise()
	// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
	// This enables us to set up storage for the data ahead of time.
	void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
		for( std::size_t i = 0; i < m_result->size(); ++i ) {
			m_result->at( i ).reserve( max_entries ) ;
		}
	}
	
	// Called once per sample to determine whether we want data for this sample
	bool set_sample( std::size_t i ) {
		m_sample_i = i ;
		// Yes, here we want info for all samples.
		return true ;
	}
	
	// Called once per sample to set the number of probabilities that are present.
	void set_number_of_entries(
		std::size_t ploidy,
		std::size_t number_of_entries,
		genfile::OrderType order_type,
		genfile::ValueType value_type
	) {
		assert( value_type == genfile::eProbability ) ;
		m_result->at( m_sample_i ).resize( number_of_entries ) ;
		m_entry_i = 0 ;
	}

	// Called once for each genotype (or haplotype) probability per sample.
	void set_value( uint32_t, double value ) {
		m_result->at( m_sample_i ).at( m_entry_i++ ) = value ;
	}

	// Ditto, but called if data is missing for this sample.
	void set_value( uint32_t, genfile::MissingValue value ) {
		// Here we encode missing probabilities with -1
		m_result->at( m_sample_i ).at( m_entry_i++ ) = -1 ;
	}

	// If present with this signature, called once after all data has been set.
	void finalise() {
		// nothing to do in this implementation.
	}

private:
	Data* m_result ;
	std::size_t m_sample_i ;
	std::size_t m_entry_i ;
} ;



// ===================================================================== //

/// Progress object
class COREARRAY_DLL_LOCAL CProgress
{
public:
	CProgress();
	CProgress(C_Int64 count, SEXP conn);

	void Reset(C_Int64 count);
	void Forward(C_Int64 val);
	void ShowProgress();

protected:
	C_Int64 fTotalCount;  ///< the total number
	C_Int64 fCounter;  ///< the current counter
	double _start, _step;
	C_Int64 _hit;
	vector< pair<double, time_t> > _timer;
	time_t _start_time, _last_time, _last_check_time;
	Rconnection progress_conn;
};

static const int PROGRESS_BAR_CHAR_NUM = 50;
static const double S_MIN  =  60;
static const double S_HOUR =  60 * S_MIN;
static const double S_DAY  =  24 * S_HOUR;
static const double S_YEAR = 365 * S_DAY;

static const char *time_str(double s)
{
	if (R_FINITE(s))
	{
		static char buffer[64];
		if (s < S_MIN)
			snprintf(buffer, sizeof(buffer), "%.0fs", s);
		else if (s < S_HOUR)
			snprintf(buffer, sizeof(buffer), "%.1fm", s/S_MIN);
		else if (s < S_DAY)
			snprintf(buffer, sizeof(buffer), "%.1fh", s/S_HOUR);
		else if (s < S_YEAR)
			snprintf(buffer, sizeof(buffer), "%.1fd", s/S_DAY);
		else
			snprintf(buffer, sizeof(buffer), "%.1f years", s/S_YEAR);
		return buffer;
	} else
		return "---";
}

extern "C"
{
	static void chkIntFn(void *dummy) { R_CheckUserInterrupt(); }
}

static bool CheckInterrupt()
{
	// this will call the above in a top-level context so it won't longjmp-out of your context
	return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

inline static void put_text(Rconnection conn, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	(*conn->vfprintf)(conn, fmt, args);
	va_end(args);
}



CProgress::CProgress()
{
	fTotalCount = 0;
	fCounter = 0;
	progress_conn = NULL;
}

CProgress::CProgress(C_Int64 count, SEXP conn)
{
	fTotalCount = 0;
	fCounter = 0;
	progress_conn = (!conn || Rf_isNull(conn)) ? NULL : R_GetConnection(conn);
	Reset(count);
}

void CProgress::Reset(C_Int64 count)
{
	bool flag = (fTotalCount==0) || (fCounter > 0);
	fTotalCount = count;
	fCounter = 0;
	if (count > 0)
	{
		int n = 100;
		if (n > count) n = count;
		if (n < 1) n = 1;
		_start = _step = (double)count / n;
		_hit = (C_Int64)(_start);
		double percent = (double)fCounter / count;
		time_t s; time(&s);
		_start_time = _last_time = _last_check_time = s;
		_timer.clear();
		_timer.reserve(128);
		_timer.push_back(pair<double, time_t>(percent, s));
		if (flag) ShowProgress();
	}
}

void CProgress::Forward(C_Int64 val)
{
	if (fTotalCount > 0)
	{
		fCounter += val;
		if (fCounter >= _hit)
		{
			do {
				_start += _step;
				_hit = (C_Int64)(_start);
			} while (fCounter >= _hit);
			if (_hit > fTotalCount) _hit = fTotalCount;
			ShowProgress();
		}
		// check whether user interrupts
		time_t now; time(&now);
		if (difftime(now, _last_check_time) > 0.25)
		{
			_last_check_time = now;
			if (CheckInterrupt())
				throw ErrCoreArray("User interrupts!");
		}
	}
}

void CProgress::ShowProgress()
{
	if (fTotalCount > 0)
	{
		char bar[PROGRESS_BAR_CHAR_NUM + 1];
		double p = (double)fCounter / fTotalCount;
		int n = (int)round(p * PROGRESS_BAR_CHAR_NUM);
		memset(bar, '.', sizeof(bar));
		memset(bar, '=', n);
		if ((fCounter > 0) && (n < PROGRESS_BAR_CHAR_NUM))
			bar[n] = '>';
		bar[PROGRESS_BAR_CHAR_NUM] = 0;

		// ETC: estimated time to complete
		n = (int)_timer.size() - 20;  // 20% as a sliding window size
		if (n < 0) n = 0;

		time_t now; time(&now);
		_timer.push_back(pair<double, time_t>(p, now));

		// in seconds
		double interval = difftime(now, _last_time);
		double s = difftime(now, _timer[n].second);
		double diff = p - _timer[n].first;
		if (diff > 0)
			s = s / diff * (1 - p);
		else
			s = R_NaN;
		p *= 100;

		// show
		_last_time = now;
		if (fCounter >= fTotalCount)
		{
			s = difftime(_last_time, _start_time);
			if (!progress_conn)
				Rprintf("\r[%s] 100%%, completed, %s\n", bar, time_str(s));
			else {
				put_text(progress_conn, "[%s] 100%%, completed, %s\n", bar, time_str(s));
				(*progress_conn->fflush)(progress_conn);
			}
		} else if ((interval >= 5) || (fCounter <= 0))
		{
			if (!progress_conn)
				Rprintf("\r[%s] %2.0f%%, ETC: %s        ", bar, p, time_str(s));
			else {
				put_text(progress_conn, "[%s] %2.0f%%, ETC: %s\n", bar, p, time_str(s));
				(*progress_conn->fflush)(progress_conn);
			}
			// fflush(stdout);
		}
	}
}



extern "C"
{

/// Return true, if matching
inline static bool StrCaseCmp(const char *prefix, const char *txt, size_t nmax)
{
	while (*prefix && *txt && nmax>0)
	{
		if (toupper(*prefix) != toupper(*txt))
			return false;
		prefix ++; txt ++; nmax --;
	}
	return (*prefix == 0);
}


/// Get the info of bgen file
COREARRAY_DLL_EXPORT SEXP SEQ_BGEN_Info(SEXP bgen_fn)
{
	// only need bgen library version information?
	if (Rf_isNull(bgen_fn)) return mkString(bgen_version);

	const char *filename = CHAR(STRING_ELT(bgen_fn, 0));
	COREARRAY_TRY

		// open the bgen file
		ifstream m_stream(filename, ifstream::binary);
		if (!m_stream)
			throw ErrCoreArray("Can't open the file '%s'.", filename);

		// read the offset, header, and sample IDs if present.
		uint32_t m_offset; // offset byte from top of bgen file
		bgen::read_offset(m_stream, &m_offset);
		// bgen::Context object holds information from the header block
		bgen::Context m_context;
		bgen::read_header_block(m_stream, &m_context);

		// output
		uint32_t v;
		const char *s;
		const size_t nlist = 9;
		PROTECT(rv_ans = NEW_LIST(nlist));
			// the number of samples
			SET_ELEMENT(rv_ans, 0,
				ScalarInteger(m_context.number_of_samples));
			// the number of variants
			SET_ELEMENT(rv_ans, 1,
				ScalarInteger(m_context.number_of_variants));
			// compression
			v = m_context.flags & bgen::e_CompressedSNPBlocks;
			if (v == bgen::e_NoCompression)
				s = "";
			else if (v == bgen::e_ZlibCompression)
				s = "zlib";
			else if (v == bgen::e_ZstdCompression)
				s = "zstd";
			else
				s = "unknown";
			SET_ELEMENT(rv_ans, 2, mkString(s));
			// layout version
			v = m_context.flags & bgen::e_Layout2;
			if (v & bgen::e_Layout2)
				s = "v1.2";
			else
				s = "v1.1";
			SET_ELEMENT(rv_ans, 3, mkString(s));

			// sample IDs
			if (m_context.flags & bgen::e_SampleIdentifiers)
			{
				vector<string> id_lst;
				bgen::read_sample_identifier_block(m_stream, m_context,
					[&](string id) { id_lst.push_back(id); });
				SEXP ss = NEW_CHARACTER(id_lst.size());
				SET_ELEMENT(rv_ans, nlist-1, ss);
				for (size_t i=0; i < id_lst.size(); i++)
					SET_STRING_ELT(ss, i, mkChar(id_lst[i].c_str()));
			}

			// pack information
			if ((m_context.flags & bgen::e_Layout2) &&
				m_context.number_of_samples>0 && m_context.number_of_variants>0)
			{
				// used for v1.2 only
				vector<genfile::byte_t> m_buffer1, m_buffer2;
				// variant info (unused)
				string chr, snpid, rsid, s;
				uint32_t position;
				bgen::read_snp_identifying_data(m_stream, m_context,
					&snpid, &rsid, &chr, &position,
					[](size_t n) { }, [](size_t i, string const& allele) { }
				);

				genfile::byte_t *buf_ptr, *buf_end;
				bgen::read_genotype_data_block(m_stream, m_context, &m_buffer1);
				if ((m_context.flags & bgen::e_CompressedSNPBlocks) !=
					bgen::e_NoCompression)
				{
					bgen::uncompress_probability_data(m_context, m_buffer1, &m_buffer2);
					buf_ptr = &m_buffer2[0];
					buf_end = &m_buffer2[0] + m_buffer2.size();
				} else {
					buf_ptr = &m_buffer1[0];
					buf_end = &m_buffer1[0] + m_buffer1.size();
				}

				bgen::v12::GenotypeDataBlock bk(m_context, buf_ptr, buf_end);
				// unphased
				SET_ELEMENT(rv_ans, 4, ScalarLogical(!bk.phased));
				// bits
				SET_ELEMENT(rv_ans, 5, ScalarInteger(bk.bits));
				// ploidy.min
				SET_ELEMENT(rv_ans, 6, ScalarInteger(bk.ploidyExtent[0]));
				// ploidy.max
				SET_ELEMENT(rv_ans, 7, ScalarInteger(bk.ploidyExtent[1]));
			} else {
				SET_ELEMENT(rv_ans, 4, ScalarLogical(TRUE)); // unphased
				SET_ELEMENT(rv_ans, 5, ScalarInteger(8));    // bits
				SET_ELEMENT(rv_ans, 6, ScalarInteger(2));    // ploidy.min
				SET_ELEMENT(rv_ans, 7, ScalarInteger(2));    // ploidy.max
			}

		UNPROTECT(1);

	COREARRAY_CATCH
}



/// Import a bgen file
COREARRAY_DLL_EXPORT SEXP SEQ_BGEN_Import(SEXP bgen_fn, SEXP gds_root, SEXP ChrPrefix,
	SEXP Start, SEXP Count, SEXP progfile, SEXP Verbose)
{
	const char *filename = CHAR(STRING_ELT(bgen_fn, 0));
	int start = Rf_asInteger(Start);
	if (start < 1) start = 1;
	int count = Rf_asInteger(Count);
	int verbose = Rf_asLogical(Verbose);
	if (verbose == NA_LOGICAL)
		Rf_error("'verbose' must be TRUE or FALSE.");

	COREARRAY_TRY

		// open the bgen file
		ifstream m_stream(filename, ifstream::binary);
		if (!m_stream)
			throw ErrCoreArray("Can't open the file '%s'.", filename);

		// read the offset, header, and sample IDs if present.
		uint32_t m_offset; // offset byte from top of bgen file
		bgen::read_offset(m_stream, &m_offset);
		// bgen::Context object holds information from the header block
		bgen::Context m_context;
		bgen::read_header_block(m_stream, &m_context);
		// jump to the first variant data block
		m_stream.seekg(m_offset + 4);

		// GDS nodes
		PdAbstractArray Root = GDS_R_SEXP2Obj(gds_root, FALSE);

		PdAbstractArray varIdx = GDS_Node_Path(Root, "variant.id", TRUE);
		PdAbstractArray varChr = GDS_Node_Path(Root, "chromosome", TRUE);
		PdAbstractArray varPos = GDS_Node_Path(Root, "position", TRUE);
		PdAbstractArray varRSID = GDS_Node_Path(Root, "annotation/id", TRUE);
		PdAbstractArray varAllele = GDS_Node_Path(Root, "allele", TRUE);

		PdAbstractArray varQual = GDS_Node_Path(Root, "annotation/qual", TRUE);
		PdAbstractArray varFilter = GDS_Node_Path(Root, "annotation/filter", TRUE);

		PdAbstractArray varGeno = GDS_Node_Path(Root, "genotype/data", FALSE);
		PdAbstractArray varGenoLen = GDS_Node_Path(Root, "genotype/@data", FALSE);
		PdAbstractArray varPhase = GDS_Node_Path(Root, "phase/data", FALSE);

		PdAbstractArray varDS = GDS_Node_Path(Root, "annotation/format/DS/data", FALSE);
		PdAbstractArray varDSLen = GDS_Node_Path(Root, "annotation/format/DS/@data", FALSE);

		PdAbstractArray varGP = GDS_Node_Path(Root, "annotation/format/GP/data", FALSE);
		PdAbstractArray varGPLen = GDS_Node_Path(Root, "annotation/format/GP/@data", FALSE);

		// initialize
		string chr, snpid, rsid, s;
		uint32_t position;
		vector<string> alleles;

		// chromosome coding prefix
		vector<const char *> ChrPref;
		for (size_t i=0; i < Rf_length(ChrPrefix); i++)
			ChrPref.push_back(CHAR(STRING_ELT(ChrPrefix, i)));

		// skip
		for (int idx=1; idx < start; idx++)
		{
			bgen::read_snp_identifying_data(m_stream, m_context,
				&snpid, &rsid, &chr, &position,
				[&alleles](size_t n) { alleles.resize(n); },
				[&alleles](size_t i, string const& allele) { alleles.at(i) = allele; }
			);
			bgen::ignore_genotype_data_block(m_stream, m_context);
		}
		// the number of variants computed
		if (count < 0)
			count = m_context.number_of_variants - start + 1;

		// buffers, these are used as working space by bgen implementation
		vector<vector<double>> probs;
		vector<genfile::byte_t> m_buffer1, m_buffer2;
		const size_t nSamp = m_context.number_of_samples;
		C_Int32 I32;
		double F64;
		vector<double> F64s(nSamp);
		vector<C_UInt8> I8s(nSamp*2);
		vector<C_UInt8> Zeros(nSamp, 0);
		// progress information
		CProgress prog((verbose || !Rf_isNull(progfile)) ? count : -1, progfile);

		// for-loop
		for (size_t idx=start; count > 0; idx++, count--)
		{
			bgen::read_snp_identifying_data(m_stream, m_context,
				&snpid, &rsid, &chr, &position,
				[&alleles](size_t n) { alleles.resize(n); },
				[&alleles](size_t i, string const& allele) { alleles.at(i) = allele; }
			);
			// variant.id
			I32 = idx;
			GDS_Array_AppendData(varIdx, 1, &I32, svInt32);

			// chromosome
			vector<const char *>::iterator ps;
			for (ps=ChrPref.begin(); ps != ChrPref.end(); ps++)
			{
				if (StrCaseCmp(*ps, chr.c_str(), chr.size()))
				{
					chr.erase(0, strlen(*ps));
					break;
				}
			}
			GDS_Array_AppendData(varChr, 1, &chr, svStrUTF8);

			// position
			I32 = position;
			GDS_Array_AppendData(varPos, 1, &I32, svInt32);
			// ID
			GDS_Array_AppendData(varRSID, 1, &rsid, svStrUTF8);
			// REF, ALT
			s = alleles[0];
			for (size_t i=1; i < alleles.size(); i++)
			{
				s.append(1, ',');
				s.append(alleles[i]);
			}
			GDS_Array_AppendData(varAllele, 1, &s, svStrUTF8);
			// QUAL
			F64 = R_NaN;
			GDS_Array_AppendData(varQual, 1, &F64, svFloat64);
			// FILTER
			I32 = NA_INTEGER;
			GDS_Array_AppendData(varFilter, 1, &I32, svInt32);

			// read probs
			ProbSetter setter(&probs);
			bgen::read_and_parse_genotype_data_block<ProbSetter>(m_stream,
				m_context, setter, &m_buffer1, &m_buffer2);

			// integer genotypes
			if (varGeno)
			{
				// write integer genotypes
				C_UInt8 *g = &I8s[0];
				for (size_t i=0; i < nSamp; i++, g+=2)
				{
					vector<double> &p = probs[i];
					if (p.size() == 3)
					{
						if (p[0]!=-1 && p[1]!=-1 && p[2]!=-1)
						{
							if (p[0] >= p[1])
							{
								g[0] = g[1] = (p[0] >= p[2]) ? 0 : 1;
							} else {
								if (p[1] > p[2])
									{ g[0] = 0; g[1] = 1; }
								else
									g[0] = g[1] = 1;
							}
						} else
							g[0] = g[1] = 0xFF;
					} else
						throw "Only support bi-allelic sites for integer genotypes.";
				}
				GDS_Array_AppendData(varGeno, nSamp*2, &I8s[0], svUInt8);
				// write data size
				I32 = 1;
				GDS_Array_AppendData(varGenoLen, 1, &I32, svInt32);
				// write phase
				GDS_Array_AppendData(varPhase, nSamp, &Zeros[0], svUInt8);
			}

			// dosages
			if (varDS)
			{
				// write dosages
				for (size_t i=0; i < nSamp; i++)
				{
					vector<double> &p = probs[i];
					if (p.size()==3 && p[0]!=-1 && p[1]!=-1 && p[2]!=-1)
						F64s[i] = 2*p[2] + p[1];
					else
						F64s[i] = R_NaN;
				}
				GDS_Array_AppendData(varDS, nSamp, &F64s[0], svFloat64);
				// write data size
				I32 = 1;
				GDS_Array_AppendData(varDSLen, 1, &I32, svInt32);
			}

			// genotype probabilities
			if (varGP)
			{
				// find the max length
				size_t max_len = 0;
				for (size_t i=0; i < nSamp; i++)
				{
					if (max_len < probs[i].size())
						max_len = probs[i].size();
				}
				// write prob
				for (size_t j=(max_len>=2)?1:0; j < max_len; j++)
				{
					for (size_t i=0; i < nSamp; i++)
					{
						vector<double> &p = probs[i];
						F64s[i] = (j<p.size() && p[j]!=-1) ? p[j] : R_NaN;
					}
					GDS_Array_AppendData(varGP, nSamp, &F64s[0], svFloat64);
				}
				// write data size
				I32 = max_len;
				GDS_Array_AppendData(varGPLen, 1, &I32, svInt32);
			}

			prog.Forward(1);
		}

	COREARRAY_CATCH
}



// ===================================================================== //
// Format conversion from GDS to BGEN (single-threaded writer)
//
// The writer keeps an output file stream open across multiple .Call()s via
// an R external pointer:
//     SEQ_BGEN_ExportBegin() -> writes the offset/header/sample-id blocks
//     SEQ_BGEN_ExportData()  -> appends a block of (bi-allelic) variants
//     SEQ_BGEN_ExportEnd()   -> flushes & closes the stream
// Only bi-allelic, diploid, unphased data (BGEN layout v1.2) is supported in
// this MVP, mirroring the bi-allelic limitation of the import path.

/// The on-disk BGEN writer state attached to an R external pointer
struct CBGenWriter
{
	std::ofstream stream;     ///< the output file stream
	bgen::Context context;    ///< the BGEN header context
	int bits;                 ///< the number of bits per probability
	bool phased;              ///< whether to write phased haplotypes (unused in MVP)
	std::string filename;     ///< the output file name
	C_Int64 nvar_written;     ///< the number of variants written so far
	std::vector<genfile::byte_t> buf_id, buf1, buf2;  ///< working buffers
};

/// Free the writer object when the external pointer is garbage-collected
static void bgen_writer_finalizer(SEXP ptr)
{
	CBGenWriter *w = (CBGenWriter *)R_ExternalPtrAddr(ptr);
	if (w)
	{
		if (w->stream.is_open()) w->stream.close();
		delete w;
		R_ClearExternalPtr(ptr);
	}
}

static CBGenWriter *get_bgen_writer(SEXP ptr)
{
	CBGenWriter *w = (CBGenWriter *)R_ExternalPtrAddr(ptr);
	if (!w) throw ErrCoreArray("The BGEN writer has been closed.");
	return w;
}


/// Create a BGEN file and write the header; return an external pointer
COREARRAY_DLL_EXPORT SEXP SEQ_BGEN_ExportBegin(SEXP OutFN, SEXP NSamp,
	SEXP NVariant, SEXP SampleID, SEXP Bits, SEXP Compress, SEXP Phased)
{
	const char *filename = CHAR(STRING_ELT(OutFN, 0));
	const int nSamp = Rf_asInteger(NSamp);
	const int nVariant = Rf_asInteger(NVariant);
	const int bits = Rf_asInteger(Bits);
	const char *compress = CHAR(STRING_ELT(Compress, 0));
	const int phased = Rf_asLogical(Phased);

	COREARRAY_TRY

		if (nSamp <= 0)
			throw ErrCoreArray("No sample to export.");
		if (bits < 1 || bits > 32)
			throw ErrCoreArray("'bits' should be between 1 and 32.");

		CBGenWriter *w = new CBGenWriter;
		w->filename = filename;
		w->bits = bits;
		w->phased = (phased == TRUE);
		w->nvar_written = 0;

		// compression flag
		uint32_t cflag;
		if (!strcmp(compress, "none") || !compress[0])
			cflag = bgen::e_NoCompression;
		else if (!strcmp(compress, "zlib"))
			cflag = bgen::e_ZlibCompression;
		else if (!strcmp(compress, "zstd"))
			cflag = bgen::e_ZstdCompression;
		else
			throw ErrCoreArray("'compress' should be 'zlib', 'zstd' or 'none'.");

		// sample identifiers (NULL for anonymized output)
		bool has_sample_id = !Rf_isNull(SampleID);
		std::vector<std::string> sample_ids;
		if (has_sample_id)
		{
			if (Rf_length(SampleID) != nSamp)
				throw ErrCoreArray("'sample.id' should match the number of samples.");
			sample_ids.reserve(nSamp);
			for (int i=0; i < nSamp; i++)
				sample_ids.push_back(CHAR(STRING_ELT(SampleID, i)));
		}

		// set up the header context (BGEN layout v1.2)
		w->context.number_of_samples = nSamp;
		w->context.number_of_variants = nVariant;
		w->context.magic = "bgen";
		w->context.free_data = "";
		w->context.flags = bgen::e_Layout2 | cflag |
			(has_sample_id ? bgen::e_SampleIdentifiers : 0);

		// compute the offset = header block + optional sample-id block
		uint32_t offset = w->context.header_size();
		if (has_sample_id)
		{
			offset += 8;
			for (int i=0; i < nSamp; i++)
				offset += 2 + sample_ids[i].size();
		}

		// open the file & write the leading blocks
		w->stream.open(filename, std::ofstream::binary);
		if (!w->stream)
			throw ErrCoreArray("Can't create the file '%s'.", filename);
		bgen::write_offset(w->stream, offset);
		bgen::write_header_block(w->stream, w->context);
		if (has_sample_id)
			bgen::write_sample_identifier_block(w->stream, w->context, sample_ids);

		// wrap in an external pointer with a finalizer
		rv_ans = R_MakeExternalPtr(w, R_NilValue, R_NilValue);
		PROTECT(rv_ans);
		R_RegisterCFinalizerEx(rv_ans, bgen_writer_finalizer, TRUE);
		UNPROTECT(1);

	COREARRAY_CATCH
}


/// Append a block of bi-allelic variants to the BGEN file
COREARRAY_DLL_EXPORT SEXP SEQ_BGEN_ExportData(SEXP Ptr, SEXP Chr, SEXP Pos,
	SEXP RSID, SEXP Ref, SEXP Alt, SEXP Prob)
{
	COREARRAY_TRY

		CBGenWriter *w = get_bgen_writer(Ptr);
		bgen::Context &ctx = w->context;
		const size_t nSamp = ctx.number_of_samples;
		const size_t nVar = Rf_xlength(Pos);

		// Prob is a numeric array, dim = c(3, nSamp, nVar), column-major
		double *prob = REAL(Prob);
		if ((size_t)Rf_xlength(Prob) != 3 * nSamp * nVar)
			throw ErrCoreArray("Internal: 'prob' has an unexpected length.");

		const uint16_t nAlleles = 2;
		for (size_t v=0; v < nVar; v++)
		{
			std::string chr = CHAR(STRING_ELT(Chr, v));
			std::string rsid = CHAR(STRING_ELT(RSID, v));
			std::string ref = CHAR(STRING_ELT(Ref, v));
			std::string alt = CHAR(STRING_ELT(Alt, v));
			uint32_t position = (uint32_t)(INTEGER(Pos)[v]);
			std::string alleles[2] = { ref, alt };

			// variant identifying data (SNPID left empty; RSID from annotation/id)
			bgen::write_snp_identifying_data(&w->buf_id, ctx,
				std::string(), rsid, chr, position, nAlleles,
				[&alleles](uint16_t i) -> std::string { return alleles[i]; });
			w->stream.write((const char *)&w->buf_id[0], w->buf_id.size());

			// genotype probability data block
			bgen::GenotypeDataBlockWriter writer(&w->buf1, &w->buf2, ctx, w->bits);
			writer.initialise(nSamp, nAlleles);
			double *p = prob + 3 * nSamp * v;
			for (size_t i=0; i < nSamp; i++, p+=3)
			{
				writer.set_sample(i);
				writer.set_number_of_entries(2, 3,
					genfile::ePerUnorderedGenotype, genfile::eProbability);
				if (p[0]!=p[0] || p[1]!=p[1] || p[2]!=p[2])  // any NaN => missing
				{
					for (uint32_t e=0; e < 3; e++)
						writer.set_value(e, genfile::MissingValue());
				} else {
					// clamp to [0,1] and renormalize so the 3 probs sum to 1
					double g[3];
					double s = 0;
					for (int e=0; e < 3; e++)
					{
						double x = p[e];
						if (x < 0) x = 0; else if (x > 1) x = 1;
						g[e] = x; s += x;
					}
					if (s > 0)
						{ g[0]/=s; g[1]/=s; g[2]/=s; }
					else
						{ g[0]=1; g[1]=g[2]=0; }
					for (uint32_t e=0; e < 3; e++)
						writer.set_value(e, g[e]);
				}
			}
			writer.finalise();
			std::pair<genfile::byte_t const*, genfile::byte_t const*> r = writer.repr();
			w->stream.write((const char *)r.first, r.second - r.first);

			w->nvar_written++;
		}

	COREARRAY_CATCH
}


/// Flush and close the BGEN file
COREARRAY_DLL_EXPORT SEXP SEQ_BGEN_ExportEnd(SEXP Ptr)
{
	COREARRAY_TRY

		CBGenWriter *w = get_bgen_writer(Ptr);
		w->stream.flush();
		if (!w->stream)
			throw ErrCoreArray("Error writing to '%s'.", w->filename.c_str());
		C_Int64 expected = w->context.number_of_variants;
		C_Int64 got = w->nvar_written;
		w->stream.close();
		// detach: free now instead of waiting for GC
		delete w;
		R_ClearExternalPtr(Ptr);
		if (got != expected)
			throw ErrCoreArray("Wrote %lld variants but the header declared %lld.",
				(long long)got, (long long)expected);

	COREARRAY_CATCH
}



/// Initialize the package
COREARRAY_DLL_EXPORT void R_init_gds2bgen(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }
	static R_CallMethodDef callMethods[] =
	{
		CALL(SEQ_BGEN_Info, 1),
		CALL(SEQ_BGEN_Import, 6),
		CALL(SEQ_BGEN_ExportBegin, 7),
		CALL(SEQ_BGEN_ExportData, 7),
		CALL(SEQ_BGEN_ExportEnd, 1),
		{ NULL, NULL, 0 }
	};
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	Init_GDS_Routines();
}

} // extern "C"
