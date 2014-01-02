static char module_doc[] =
"This module provides a BLAS optimized\nmatrix multiply, inner product and dot for numpy arrays";

#include "Python.h"
#include "numpy/noprefix.h"
#ifndef CBLAS_HEADER
#define CBLAS_HEADER "cblas.h"
#endif
#include CBLAS_HEADER

#include <stdio.h>

static PyArray_DotFunc *oldFunctions[PyArray_NTYPES];

static void
FLOAT_dot(void *a, intp stridea, void *b, intp strideb, void *res,
	  intp n, void *tmp)
{
    register int na = stridea / sizeof(float);
    register int nb = strideb / sizeof(float);

    if ((sizeof(float) * na == stridea) &&
	(sizeof(float) * nb == strideb) &&
	(na >= 0) && (nb >= 0))
	    *((float *)res) = cblas_sdot((int)n, (float *)a, na, (float *)b, nb);

    else
	    oldFunctions[PyArray_FLOAT](a, stridea, b, strideb, res, n, tmp);
}

static void
DOUBLE_dot(void *a, intp stridea, void *b, intp strideb, void *res,
	   intp n, void *tmp)
{
    register int na = stridea / sizeof(double);
    register int nb = strideb / sizeof(double);

    if ((sizeof(double) * na == stridea) &&
	(sizeof(double) * nb == strideb) &&
	(na >= 0) && (nb >= 0))
	    *((double *)res) = cblas_ddot((int)n, (double *)a, na, (double *)b, nb);
    else
	    oldFunctions[PyArray_DOUBLE](a, stridea, b, strideb, res, n, tmp);
}

static void
CFLOAT_dot(void *a, intp stridea, void *b, intp strideb, void *res,
	   intp n, void *tmp)
{

    register int na = stridea / sizeof(cfloat);
    register int nb = strideb / sizeof(cfloat);

    if ((sizeof(cfloat) * na == stridea) &&
	(sizeof(cfloat) * nb == strideb) &&
	(na >= 0) && (nb >= 0))
	    cblas_cdotu_sub((int)n, (float *)a, na, (float *)b, nb, (float *)res);
    else
	    oldFunctions[PyArray_CFLOAT](a, stridea, b, strideb, res, n, tmp);
}

static void
CDOUBLE_dot(void *a, intp stridea, void *b, intp strideb, void *res,
	    intp n, void *tmp)
{
    register int na = stridea / sizeof(cdouble);
    register int nb = strideb / sizeof(cdouble);

    if ((sizeof(cdouble) * na == stridea) &&
	(sizeof(cdouble) * nb == strideb) &&
	(na >= 0) && (nb >= 0))
	    cblas_zdotu_sub((int)n, (double *)a, na, (double *)b, nb, (double *)res);
    else
	    oldFunctions[PyArray_CDOUBLE](a, stridea, b, strideb, res, n, tmp);
}


static Bool altered=FALSE;

static char doc_alterdot[] = "alterdot() changes all dot functions to use blas.";

static PyObject *
dotblas_alterdot(PyObject *dummy, PyObject *args)
{
    PyArray_Descr *descr;

    if (!PyArg_ParseTuple(args, "")) return NULL;

    /* Replace the dot functions to the ones using blas */

    if (!altered) {
	descr = PyArray_DescrFromType(PyArray_FLOAT);
	oldFunctions[PyArray_FLOAT] = descr->f->dotfunc;
	descr->f->dotfunc = (PyArray_DotFunc *)FLOAT_dot;

	descr = PyArray_DescrFromType(PyArray_DOUBLE);
	oldFunctions[PyArray_DOUBLE] = descr->f->dotfunc;
	descr->f->dotfunc = (PyArray_DotFunc *)DOUBLE_dot;

	descr = PyArray_DescrFromType(PyArray_CFLOAT);
	oldFunctions[PyArray_CFLOAT] = descr->f->dotfunc;
	descr->f->dotfunc = (PyArray_DotFunc *)CFLOAT_dot;

	descr = PyArray_DescrFromType(PyArray_CDOUBLE);
	oldFunctions[PyArray_CDOUBLE] = descr->f->dotfunc;
	descr->f->dotfunc = (PyArray_DotFunc *)CDOUBLE_dot;

	altered = TRUE;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static char doc_restoredot[] = "restoredot() restores dots to defaults.";

static PyObject *
dotblas_restoredot(PyObject *dummy, PyObject *args)
{
    PyArray_Descr *descr;

    if (!PyArg_ParseTuple(args, "")) return NULL;

    if (altered) {
	descr = PyArray_DescrFromType(PyArray_FLOAT);
	descr->f->dotfunc = oldFunctions[PyArray_FLOAT];
	oldFunctions[PyArray_FLOAT] = NULL;
	Py_XDECREF(descr);

	descr = PyArray_DescrFromType(PyArray_DOUBLE);
	descr->f->dotfunc = oldFunctions[PyArray_DOUBLE];
	oldFunctions[PyArray_DOUBLE] = NULL;
	Py_XDECREF(descr);

	descr = PyArray_DescrFromType(PyArray_CFLOAT);
	descr->f->dotfunc = oldFunctions[PyArray_CFLOAT];
	oldFunctions[PyArray_CFLOAT] = NULL;
	Py_XDECREF(descr);

	descr = PyArray_DescrFromType(PyArray_CDOUBLE);
	descr->f->dotfunc = oldFunctions[PyArray_CDOUBLE];
	oldFunctions[PyArray_CDOUBLE] = NULL;
	Py_XDECREF(descr);

	altered = FALSE;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

typedef enum {_scalar, _column, _row, _matrix} MatrixShape;

static MatrixShape
_select_matrix_shape(PyArrayObject *array)
{
    switch (array->nd) {
    case 0:
	return _scalar;
    case 1:
	if (array->dimensions[0] > 1)
	    return _column;
	return _scalar;
    case 2:
	if (array->dimensions[0] > 1) {
	    if (array->dimensions[1] == 1)
		return _column;
	    else
		return _matrix;
	}
	if (array->dimensions[1] == 1)
	    return _scalar;
	return _row;
    }
    return _matrix;
}


/* This also makes sure that the data segment is aligned with
   an itemsize address as well by returning one if not true. 
*/
static int
_bad_strides(PyArrayObject *ap)
{
    register int itemsize = PyArray_ITEMSIZE(ap);
    register int i, N=PyArray_NDIM(ap);
    register intp *strides = PyArray_STRIDES(ap);

    if (((intp)(ap->data) % itemsize) != 0)
	return 1;
    for (i=0; i<N; i++) {
	if ((strides[i] < 0) || (strides[i] % itemsize) != 0) 
	    return 1;
    }

    return 0;  
}

static char doc_matrixproduct[] = "dot(a,b)\nReturns the dot product of a and b for arrays of floating point types.\nLike the generic numpy equivalent the product sum is over\nthe last dimension of a and the second-to-last dimension of b.\nNB: The first argument is not conjugated.";

static PyObject *
dotblas_matrixproduct(PyObject *dummy, PyObject *args)
{
    PyObject *op1, *op2;
    PyArrayObject *ap1=NULL, *ap2=NULL, *ret=NULL;
    int j, l, lda, ldb, ldc;
    int typenum, nd;
    intp ap1stride=0;
    intp dimensions[MAX_DIMS];
    intp numbytes;
    static const float oneF[2] = {1.0, 0.0};
    static const float zeroF[2] = {0.0, 0.0};
    static const double oneD[2] = {1.0, 0.0};
    static const double zeroD[2] = {0.0, 0.0};
    double prior1, prior2;
    PyTypeObject *subtype;
    PyArray_Descr *dtype;
    MatrixShape ap1shape, ap2shape;

    if (!PyArg_ParseTuple(args, "OO", &op1, &op2)) return NULL;

    /*
     * "Matrix product" using the BLAS.
     * Only works for float double and complex types.
     */

    typenum = PyArray_ObjectType(op1, 0);
    typenum = PyArray_ObjectType(op2, typenum);

    /* This function doesn't handle other types */
    if ((typenum != PyArray_DOUBLE && typenum != PyArray_CDOUBLE &&
	 typenum != PyArray_FLOAT && typenum != PyArray_CFLOAT)) {
	return PyArray_Return((PyArrayObject *)PyArray_MatrixProduct(op1, op2));
    }

    dtype = PyArray_DescrFromType(typenum);
    ap1 = (PyArrayObject *)PyArray_FromAny(op1, dtype, 0, 0, ALIGNED, NULL);
    if (ap1 == NULL) return NULL;
    Py_INCREF(dtype);
    ap2 = (PyArrayObject *)PyArray_FromAny(op2, dtype, 0, 0, ALIGNED, NULL);
    if (ap2 == NULL) goto fail;


    if ((ap1->nd > 2) || (ap2->nd > 2)) {
	/* This function doesn't handle dimensions greater than 2 
	   (or negative striding)  -- other
	   than to ensure the dot function is altered
	*/
	if (!altered) {
	    /* need to alter dot product */
	    PyObject *tmp1, *tmp2;
	    tmp1 = PyTuple_New(0);
	    tmp2 = dotblas_alterdot(NULL, tmp1);
	    Py_DECREF(tmp1);
	    Py_DECREF(tmp2);
	}
	ret = (PyArrayObject *)PyArray_MatrixProduct((PyObject *)ap1,
						     (PyObject *)ap2);
	Py_DECREF(ap1);
	Py_DECREF(ap2);
	return PyArray_Return(ret);
    }

    if (_bad_strides(ap1)) {
	    op1 = PyArray_NewCopy(ap1, PyArray_ANYORDER);
	    Py_DECREF(ap1);
	    ap1 = (PyArrayObject *)op1;
	    if (ap1 == NULL) goto fail;
    }
    if (_bad_strides(ap2)) {
	    op2 = PyArray_NewCopy(ap2, PyArray_ANYORDER);
	    Py_DECREF(ap2);
	    ap2 = (PyArrayObject *)op2;
	    if (ap2 == NULL) goto fail;
    }
    ap1shape = _select_matrix_shape(ap1);
    ap2shape = _select_matrix_shape(ap2);

    if (ap1shape == _scalar || ap2shape == _scalar) {
        PyArrayObject *oap1, *oap2;
        oap1 = ap1; oap2 = ap2;
	/* One of ap1 or ap2 is a scalar */
	if (ap1shape == _scalar) { 		/* Make ap2 the scalar */
	    PyArrayObject *t = ap1;
	    ap1 = ap2;
	    ap2 = t;
	    ap1shape = ap2shape;
	    ap2shape = _scalar;
	}

	if (ap1shape == _row) ap1stride = ap1->strides[1];
	else if (ap1->nd > 0) ap1stride = ap1->strides[0];

 	if (ap1->nd == 0 || ap2->nd == 0) {
            intp *thisdims;
            if (ap1->nd == 0) {
                nd = ap2->nd;
                thisdims = ap2->dimensions;
            }
            else {
                nd = ap1->nd;
                thisdims = ap1->dimensions;
            }
            l = 1;
            for (j=0; j<nd; j++) {
                dimensions[j] = thisdims[j];
                l *= dimensions[j];
            }
        }
        else {
            l = oap1->dimensions[oap1->nd-1];

            if (oap2->dimensions[0] != l) {
                PyErr_SetString(PyExc_ValueError, "matrices are not aligned");
                goto fail;
            }
            nd = ap1->nd + ap2->nd - 2;
            /* nd = 0 or 1 or 2 */
            /* If nd == 0 do nothing ... */
            if (nd == 1) {
                /* Either ap1->nd is 1 dim or ap2->nd is 1 dim
                   and the other is 2-dim */
                dimensions[0] = (oap1->nd == 2) ? oap1->dimensions[0] : oap2->dimensions[1];
                l = dimensions[0];
                /* Fix it so that dot(shape=(N,1), shape=(1,))
                   and dot(shape=(1,), shape=(1,N)) both return
                   an (N,) array (but use the fast scalar code)
                */
            }
            else if (nd == 2) {
                dimensions[0] = oap1->dimensions[0];
                dimensions[1] = oap2->dimensions[1];
		/* We need to make sure that dot(shape=(1,1), shape=(1,N))
		   and dot(shape=(N,1),shape=(1,1)) uses
		   scalar multiplication appropriately
		*/
		if (ap1shape == _row) l = dimensions[1];
		else l = dimensions[0];
            }
	}
    }
    else { /* (ap1->nd <= 2 && ap2->nd <= 2) */
	/*  Both ap1 and ap2 are vectors or matrices */
	l = ap1->dimensions[ap1->nd-1];

	if (ap2->dimensions[0] != l) {
	    PyErr_SetString(PyExc_ValueError, "matrices are not aligned");
	    goto fail;
	}
	nd = ap1->nd+ap2->nd-2;

	if (nd == 1)
	    dimensions[0] = (ap1->nd == 2) ? ap1->dimensions[0] : ap2->dimensions[1];
	else if (nd == 2) {
	    dimensions[0] = ap1->dimensions[0];
	    dimensions[1] = ap2->dimensions[1];
	}
    }

    /* Choose which subtype to return */
    if (ap1->ob_type != ap2->ob_type) {
	prior2 = PyArray_GetPriority((PyObject *)ap2, 0.0);
	prior1 = PyArray_GetPriority((PyObject *)ap1, 0.0);
	subtype = (prior2 > prior1 ? ap2->ob_type : ap1->ob_type);
    }
    else {
	prior1 = prior2 = 0.0;
	subtype = ap1->ob_type;
    }

    ret = (PyArrayObject *)PyArray_New(subtype, nd, dimensions,
				       typenum, NULL, NULL, 0, 0,
				       (PyObject *)
				       (prior2 > prior1 ? ap2 : ap1));

    if (ret == NULL) goto fail;
    numbytes = PyArray_NBYTES(ret);
    memset(ret->data, 0, numbytes);
    if (numbytes==0 || l == 0) {
	    Py_DECREF(ap1);
	    Py_DECREF(ap2);
	    return PyArray_Return(ret);
    }


    if (ap2shape == _scalar) {
	/* Multiplication by a scalar -- Level 1 BLAS */
	/* if ap1shape is a matrix and we are not contiguous, then we can't
	   just blast through the entire array using a single
	   striding factor */
	NPY_BEGIN_ALLOW_THREADS

	if (typenum == PyArray_DOUBLE) {
	    if (l == 1) {
		*((double *)ret->data) = *((double *)ap2->data) * \
		    *((double *)ap1->data);
	    }
	    else if (ap1shape != _matrix) {
		cblas_daxpy(l, *((double *)ap2->data), (double *)ap1->data,
			    ap1stride/sizeof(double), (double *)ret->data, 1);
	    }
	    else {
		int maxind, oind, i, a1s, rets;
		char *ptr, *rptr;
		double val;
		maxind = (ap1->dimensions[0] >= ap1->dimensions[1] ? 0 : 1);
		oind = 1-maxind;
		ptr = ap1->data;
		rptr = ret->data;
		l = ap1->dimensions[maxind];
		val = *((double *)ap2->data);
		a1s = ap1->strides[maxind] / sizeof(double);
		rets = ret->strides[maxind] / sizeof(double);
		for (i=0; i < ap1->dimensions[oind]; i++) {
		    cblas_daxpy(l, val, (double *)ptr, a1s,
				(double *)rptr, rets);
		    ptr += ap1->strides[oind];
		    rptr += ret->strides[oind];
		}
	    }
	}
	else if (typenum == PyArray_CDOUBLE) {
	    if (l == 1) {
		cdouble *ptr1, *ptr2, *res;
		ptr1 = (cdouble *)ap2->data;
		ptr2 = (cdouble *)ap1->data;
		res = (cdouble *)ret->data;
		res->real = ptr1->real * ptr2->real - ptr1->imag * ptr2->imag;
		res->imag = ptr1->real * ptr2->imag + ptr1->imag * ptr2->real;
	    }
	    else if (ap1shape != _matrix) {
		cblas_zaxpy(l, (double *)ap2->data, (double *)ap1->data,
			    ap1stride/sizeof(cdouble), (double *)ret->data, 1);
	    }
	    else {
		int maxind, oind, i, a1s, rets;
		char *ptr, *rptr;
		double *pval;
		maxind = (ap1->dimensions[0] >= ap1->dimensions[1] ? 0 : 1);
		oind = 1-maxind;
		ptr = ap1->data;
		rptr = ret->data;
		l = ap1->dimensions[maxind];
		pval = (double *)ap2->data;
		a1s = ap1->strides[maxind] / sizeof(cdouble);
		rets = ret->strides[maxind] / sizeof(cdouble);
		for (i=0; i < ap1->dimensions[oind]; i++) {
		    cblas_zaxpy(l, pval, (double *)ptr, a1s,
				(double *)rptr, rets);
		    ptr += ap1->strides[oind];
		    rptr += ret->strides[oind];
		}
	    }
	}
	else if (typenum == PyArray_FLOAT) {
	    if (l == 1) {
		*((float *)ret->data) = *((float *)ap2->data) * \
		    *((float *)ap1->data);
	    }
	    else if (ap1shape != _matrix) {
		cblas_saxpy(l, *((float *)ap2->data), (float *)ap1->data,
			    ap1stride/sizeof(float), (float *)ret->data, 1);
	    }
	    else {
		int maxind, oind, i, a1s, rets;
		char *ptr, *rptr;
		float val;
		maxind = (ap1->dimensions[0] >= ap1->dimensions[1] ? 0 : 1);
		oind = 1-maxind;
		ptr = ap1->data;
		rptr = ret->data;
		l = ap1->dimensions[maxind];
		val = *((float *)ap2->data);
		a1s = ap1->strides[maxind] / sizeof(float);
		rets = ret->strides[maxind] / sizeof(float);
		for (i=0; i < ap1->dimensions[oind]; i++) {
		    cblas_saxpy(l, val, (float *)ptr, a1s,
				(float *)rptr, rets);
		    ptr += ap1->strides[oind];
		    rptr += ret->strides[oind];
		}
	    }
	}
	else if (typenum == PyArray_CFLOAT) {
	    if (l == 1) {
		cfloat *ptr1, *ptr2, *res;
		ptr1 = (cfloat *)ap2->data;
		ptr2 = (cfloat *)ap1->data;
		res = (cfloat *)ret->data;
		res->real = ptr1->real * ptr2->real - ptr1->imag * ptr2->imag;
		res->imag = ptr1->real * ptr2->imag + ptr1->imag * ptr2->real;
	    }
	    else if (ap1shape != _matrix) {
		cblas_caxpy(l, (float *)ap2->data, (float *)ap1->data,
			    ap1stride/sizeof(cfloat), (float *)ret->data, 1);
	    }
	    else {
		int maxind, oind, i, a1s, rets;
		char *ptr, *rptr;
		float *pval;
		maxind = (ap1->dimensions[0] >= ap1->dimensions[1] ? 0 : 1);
		oind = 1-maxind;
		ptr = ap1->data;
		rptr = ret->data;
		l = ap1->dimensions[maxind];
		pval = (float *)ap2->data;
		a1s = ap1->strides[maxind] / sizeof(cfloat);
		rets = ret->strides[maxind] / sizeof(cfloat);
		for (i=0; i < ap1->dimensions[oind]; i++) {
		    cblas_caxpy(l, pval, (float *)ptr, a1s,
				(float *)rptr, rets);
		    ptr += ap1->strides[oind];
		    rptr += ret->strides[oind];
		}
	    }
	}
	NPY_END_ALLOW_THREADS
    }
    else if ((ap2shape == _column) && (ap1shape != _matrix)) {
	int ap1s, ap2s;
	NPY_BEGIN_ALLOW_THREADS

	ap2s = ap2->strides[0] / ap2->descr->elsize;
	if (ap1shape == _row) {
	    ap1s = ap1->strides[1] / ap1->descr->elsize;
	}
	else {
	    ap1s = ap1->strides[0] / ap1->descr->elsize;
	}

	/* Dot product between two vectors -- Level 1 BLAS */
	if (typenum == PyArray_DOUBLE) {
	    double result = cblas_ddot(l, (double *)ap1->data, ap1s,
				       (double *)ap2->data, ap2s);
	    *((double *)ret->data) = result;
	}
	else if (typenum == PyArray_FLOAT) {
	    float result = cblas_sdot(l, (float *)ap1->data, ap1s,
				      (float *)ap2->data, ap2s);
	    *((float *)ret->data) = result;
	}
	else if (typenum == PyArray_CDOUBLE) {
	    cblas_zdotu_sub(l, (double *)ap1->data, ap1s,
			    (double *)ap2->data, ap2s, (double *)ret->data);
	}
	else if (typenum == PyArray_CFLOAT) {
	    cblas_cdotu_sub(l, (float *)ap1->data, ap1s,
			    (float *)ap2->data, ap2s, (float *)ret->data);
	}
	NPY_END_ALLOW_THREADS
    }
    else if (ap1shape == _matrix && ap2shape != _matrix) {
	/* Matrix vector multiplication -- Level 2 BLAS */
	/* lda must be MAX(M,1) */
	enum CBLAS_ORDER Order;
	int ap2s;

	if (!PyArray_ISONESEGMENT(ap1)) {
	    PyObject *new;
	    new = PyArray_Copy(ap1);
	    Py_DECREF(ap1);
	    ap1 = (PyArrayObject *)new;
	    if (new == NULL) goto fail;
	}
	NPY_BEGIN_ALLOW_THREADS
	if (PyArray_ISCONTIGUOUS(ap1)) {
            Order = CblasRowMajor;
            lda = (ap1->dimensions[1] > 1 ? ap1->dimensions[1] : 1);
	}
	else {
	    Order = CblasColMajor;
            lda = (ap1->dimensions[0] > 1 ? ap1->dimensions[0] : 1);
	}
	ap2s = ap2->strides[0] / ap2->descr->elsize;
	if (typenum == PyArray_DOUBLE) {
	    cblas_dgemv(Order, CblasNoTrans,
			ap1->dimensions[0], ap1->dimensions[1],
			1.0, (double *)ap1->data, lda,
			(double *)ap2->data, ap2s, 0.0, (double *)ret->data, 1);
	}
	else if (typenum == PyArray_FLOAT) {
	    cblas_sgemv(Order, CblasNoTrans,
			ap1->dimensions[0], ap1->dimensions[1],
			1.0, (float *)ap1->data, lda,
			(float *)ap2->data, ap2s, 0.0, (float *)ret->data, 1);
	}
	else if (typenum == PyArray_CDOUBLE) {
	    cblas_zgemv(Order,
			CblasNoTrans,  ap1->dimensions[0], ap1->dimensions[1],
			oneD, (double *)ap1->data, lda,
			(double *)ap2->data, ap2s, zeroD,
			(double *)ret->data, 1);
	}
	else if (typenum == PyArray_CFLOAT) {
	    cblas_cgemv(Order,
			CblasNoTrans,  ap1->dimensions[0], ap1->dimensions[1],
			oneF, (float *)ap1->data, lda,
			(float *)ap2->data, ap2s, zeroF,
			(float *)ret->data, 1);
	}
	NPY_END_ALLOW_THREADS
    }
    else if (ap1shape != _matrix && ap2shape == _matrix) {
	/* Vector matrix multiplication -- Level 2 BLAS */
	enum CBLAS_ORDER Order;
	int ap1s;

	if (!PyArray_ISONESEGMENT(ap2)) {
	    PyObject *new;
	    new = PyArray_Copy(ap2);
	    Py_DECREF(ap2);
	    ap2 = (PyArrayObject *)new;
	    if (new == NULL) goto fail;
	}
	NPY_BEGIN_ALLOW_THREADS
	if (PyArray_ISCONTIGUOUS(ap2)) {
	    Order = CblasRowMajor;
            lda = (ap2->dimensions[1] > 1 ? ap2->dimensions[1] : 1);
	}
	else {
	    Order = CblasColMajor;
            lda = (ap2->dimensions[0] > 1 ? ap2->dimensions[0] : 1);
	}
	if (ap1shape == _row) {
	    ap1s = ap1->strides[1] / ap1->descr->elsize;
	}
	else {
	    ap1s = ap1->strides[0] / ap1->descr->elsize;
	}
	if (typenum == PyArray_DOUBLE) {
	    cblas_dgemv(Order,
			CblasTrans,  ap2->dimensions[0], ap2->dimensions[1],
			1.0, (double *)ap2->data, lda,
			(double *)ap1->data, ap1s, 0.0, (double *)ret->data, 1);
	}
	else if (typenum == PyArray_FLOAT) {
	    cblas_sgemv(Order,
			CblasTrans,  ap2->dimensions[0], ap2->dimensions[1],
			1.0, (float *)ap2->data, lda,
			(float *)ap1->data, ap1s, 0.0, (float *)ret->data, 1);
	}
	else if (typenum == PyArray_CDOUBLE) {
	    cblas_zgemv(Order,
			CblasTrans,  ap2->dimensions[0], ap2->dimensions[1],
			oneD, (double *)ap2->data, lda,
			(double *)ap1->data, ap1s, zeroD, (double *)ret->data, 1);
	}
	else if (typenum == PyArray_CFLOAT) {
	    cblas_cgemv(Order,
			CblasTrans,  ap2->dimensions[0], ap2->dimensions[1],
			oneF, (float *)ap2->data, lda,
			(float *)ap1->data, ap1s, zeroF, (float *)ret->data, 1);
	}
	NPY_END_ALLOW_THREADS
    }
    else { /* (ap1->nd == 2 && ap2->nd == 2) */
	/* Matrix matrix multiplication -- Level 3 BLAS */
	/*  L x M  multiplied by M x N */
	enum CBLAS_ORDER Order;
	enum CBLAS_TRANSPOSE Trans1, Trans2;
	int M, N, L;

	/* Optimization possible: */
	/* We may be able to handle single-segment arrays here
	   using appropriate values of Order, Trans1, and Trans2.
	*/

 	if (!PyArray_ISCONTIGUOUS(ap2)) {
	    PyObject *new;
	    new = PyArray_Copy(ap2);
	    Py_DECREF(ap2);
	    ap2 = (PyArrayObject *)new;
	    if (new == NULL) goto fail;
	}
	if (!PyArray_ISCONTIGUOUS(ap1)) {
	    PyObject *new;
	    new = PyArray_Copy(ap1);
	    Py_DECREF(ap1);
	    ap1 = (PyArrayObject *)new;
	    if (new == NULL) goto fail;
	}

	NPY_BEGIN_ALLOW_THREADS

	Order = CblasRowMajor;
	Trans1 = CblasNoTrans;
	Trans2 = CblasNoTrans;
	L = ap1->dimensions[0];
	N = ap2->dimensions[1];
	M = ap2->dimensions[0];
	lda = (ap1->dimensions[1] > 1 ? ap1->dimensions[1] : 1);
	ldb = (ap2->dimensions[1] > 1 ? ap2->dimensions[1] : 1);
	ldc = (ret->dimensions[1] > 1 ? ret->dimensions[1] : 1);
	if (typenum == PyArray_DOUBLE) {
	    cblas_dgemm(Order, Trans1, Trans2,
			L, N, M,
			1.0, (double *)ap1->data, lda,
			(double *)ap2->data, ldb,
			0.0, (double *)ret->data, ldc);
	}
	else if (typenum == PyArray_FLOAT) {
	    cblas_sgemm(Order, Trans1, Trans2,
			L, N, M,
			1.0, (float *)ap1->data, lda,
			(float *)ap2->data, ldb,
			0.0, (float *)ret->data, ldc);
	}
	else if (typenum == PyArray_CDOUBLE) {
	    cblas_zgemm(Order, Trans1, Trans2,
			L, N, M,
			oneD, (double *)ap1->data, lda,
			(double *)ap2->data, ldb,
			zeroD, (double *)ret->data, ldc);
	}
	else if (typenum == PyArray_CFLOAT) {
	    cblas_cgemm(Order, Trans1, Trans2,
			L, N, M,
			oneF, (float *)ap1->data, lda,
			(float *)ap2->data, ldb,
			zeroF, (float *)ret->data, ldc);
	}
	NPY_END_ALLOW_THREADS
    }


    Py_DECREF(ap1);
    Py_DECREF(ap2);
    return PyArray_Return(ret);

 fail:
    Py_XDECREF(ap1);
    Py_XDECREF(ap2);
    Py_XDECREF(ret);
    return NULL;
}


static char doc_innerproduct[] = "innerproduct(a,b)\nReturns the inner product of a and b for arrays of floating point types.\nLike the generic NumPy equivalent the product sum is over\nthe last dimension of a and b.\nNB: The first argument is not conjugated.";

static PyObject *
dotblas_innerproduct(PyObject *dummy, PyObject *args)
{
    PyObject *op1, *op2;
    PyArrayObject *ap1, *ap2, *ret;
    int j, l, lda, ldb, ldc;
    int typenum, nd;
    intp dimensions[MAX_DIMS];
    static const float oneF[2] = {1.0, 0.0};
    static const float zeroF[2] = {0.0, 0.0};
    static const double oneD[2] = {1.0, 0.0};
    static const double zeroD[2] = {0.0, 0.0};
    PyTypeObject *subtype;
    double prior1, prior2;

    if (!PyArg_ParseTuple(args, "OO", &op1, &op2)) return NULL;

    /*
     * Inner product using the BLAS.  The product sum is taken along the last
     * dimensions of the two arrays.
     * Only speeds things up for float double and complex types.
     */


    typenum = PyArray_ObjectType(op1, 0);
    typenum = PyArray_ObjectType(op2, typenum);

    /* This function doesn't handle other types */
    if ((typenum != PyArray_DOUBLE && typenum != PyArray_CDOUBLE &&
	 typenum != PyArray_FLOAT && typenum != PyArray_CFLOAT)) {
            return PyArray_Return((PyArrayObject *)PyArray_InnerProduct(op1, op2));
    }

    ret = NULL;
    ap1 = (PyArrayObject *)PyArray_ContiguousFromObject(op1, typenum, 0, 0);
    if (ap1 == NULL) return NULL;
    ap2 = (PyArrayObject *)PyArray_ContiguousFromObject(op2, typenum, 0, 0);
    if (ap2 == NULL) goto fail;

    if ((ap1->nd > 2) || (ap2->nd > 2)) {
	/* This function doesn't handle dimensions greater than 2 -- other
	   than to ensure the dot function is altered
	*/
	if (!altered) {
	    /* need to alter dot product */
	    PyObject *tmp1, *tmp2;
	    tmp1 = PyTuple_New(0);
	    tmp2 = dotblas_alterdot(NULL, tmp1);
	    Py_DECREF(tmp1);
	    Py_DECREF(tmp2);
	}
	ret = (PyArrayObject *)PyArray_InnerProduct((PyObject *)ap1,
						    (PyObject *)ap2);
	Py_DECREF(ap1);
	Py_DECREF(ap2);
	return PyArray_Return(ret);
    }

    if (ap1->nd == 0 || ap2->nd == 0) {
	/* One of ap1 or ap2 is a scalar */
	if (ap1->nd == 0) {		/* Make ap2 the scalar */
	    PyArrayObject *t = ap1;
	    ap1 = ap2;
	    ap2 = t;
	}
	for (l = 1, j = 0; j < ap1->nd; j++) {
	    dimensions[j] = ap1->dimensions[j];
	    l *= dimensions[j];
	}
	nd = ap1->nd;
    }
    else { /* (ap1->nd <= 2 && ap2->nd <= 2) */
	/*  Both ap1 and ap2 are vectors or matrices */
	l = ap1->dimensions[ap1->nd-1];

	if (ap2->dimensions[ap2->nd-1] != l) {
	    PyErr_SetString(PyExc_ValueError, "matrices are not aligned");
	    goto fail;
	}
	nd = ap1->nd+ap2->nd-2;

	if (nd == 1)
	    dimensions[0] = (ap1->nd == 2) ? ap1->dimensions[0] : ap2->dimensions[0];
	else if (nd == 2) {
	    dimensions[0] = ap1->dimensions[0];
	    dimensions[1] = ap2->dimensions[0];
	}
    }

    /* Choose which subtype to return */
    prior2 = PyArray_GetPriority((PyObject *)ap2, 0.0);
    prior1 = PyArray_GetPriority((PyObject *)ap1, 0.0);
    subtype = (prior2 > prior1 ? ap2->ob_type : ap1->ob_type);

    ret = (PyArrayObject *)PyArray_New(subtype, nd, dimensions,
				       typenum, NULL, NULL, 0, 0,
				       (PyObject *)\
				       (prior2 > prior1 ? ap2 : ap1));

    if (ret == NULL) goto fail;
    NPY_BEGIN_ALLOW_THREADS
    memset(ret->data, 0, PyArray_NBYTES(ret));

    if (ap2->nd == 0) {
	/* Multiplication by a scalar -- Level 1 BLAS */
	if (typenum == PyArray_DOUBLE) {
	    cblas_daxpy(l, *((double *)ap2->data), (double *)ap1->data, 1,
			(double *)ret->data, 1);
	}
	else if (typenum == PyArray_CDOUBLE) {
	    cblas_zaxpy(l, (double *)ap2->data, (double *)ap1->data, 1,
			(double *)ret->data, 1);
	}
	else if (typenum == PyArray_FLOAT) {
	    cblas_saxpy(l, *((float *)ap2->data), (float *)ap1->data, 1,
			(float *)ret->data, 1);
	}
	else if (typenum == PyArray_CFLOAT) {
	    cblas_caxpy(l, (float *)ap2->data, (float *)ap1->data, 1,
			(float *)ret->data, 1);
	}
    }
    else if (ap1->nd == 1 && ap2->nd == 1) {
	/* Dot product between two vectors -- Level 1 BLAS */
	if (typenum == PyArray_DOUBLE) {
	    double result = cblas_ddot(l, (double *)ap1->data, 1,
				       (double *)ap2->data, 1);
	    *((double *)ret->data) = result;
	}
	else if (typenum == PyArray_CDOUBLE) {
	    cblas_zdotu_sub(l, (double *)ap1->data, 1,
			    (double *)ap2->data, 1, (double *)ret->data);
	}
	else if (typenum == PyArray_FLOAT) {
	    float result = cblas_sdot(l, (float *)ap1->data, 1,
				      (float *)ap2->data, 1);
	    *((float *)ret->data) = result;
	}
	else if (typenum == PyArray_CFLOAT) {
	    cblas_cdotu_sub(l, (float *)ap1->data, 1,
			    (float *)ap2->data, 1, (float *)ret->data);
	}
    }
    else if (ap1->nd == 2 && ap2->nd == 1) {
	/* Matrix-vector multiplication -- Level 2 BLAS */
	lda = (ap1->dimensions[1] > 1 ? ap1->dimensions[1] : 1);
	if (typenum == PyArray_DOUBLE) {
	    cblas_dgemv(CblasRowMajor,
			CblasNoTrans,  ap1->dimensions[0], ap1->dimensions[1],
			1.0, (double *)ap1->data, lda,
			(double *)ap2->data, 1, 0.0, (double *)ret->data, 1);
	}
	else if (typenum == PyArray_CDOUBLE) {
	    cblas_zgemv(CblasRowMajor,
			CblasNoTrans,  ap1->dimensions[0], ap1->dimensions[1],
			oneD, (double *)ap1->data, lda,
			(double *)ap2->data, 1, zeroD, (double *)ret->data, 1);
	}
	else if (typenum == PyArray_FLOAT) {
	    cblas_sgemv(CblasRowMajor,
			CblasNoTrans,  ap1->dimensions[0], ap1->dimensions[1],
			1.0, (float *)ap1->data, lda,
			(float *)ap2->data, 1, 0.0, (float *)ret->data, 1);
	}
	else if (typenum == PyArray_CFLOAT) {
	    cblas_cgemv(CblasRowMajor,
			CblasNoTrans,  ap1->dimensions[0], ap1->dimensions[1],
			oneF, (float *)ap1->data, lda,
			(float *)ap2->data, 1, zeroF, (float *)ret->data, 1);
	}
    }
    else if (ap1->nd == 1 && ap2->nd == 2) {
	/* Vector matrix multiplication -- Level 2 BLAS */
	lda = (ap2->dimensions[1] > 1 ? ap2->dimensions[1] : 1);
	if (typenum == PyArray_DOUBLE) {
	    cblas_dgemv(CblasRowMajor,
			CblasNoTrans,  ap2->dimensions[0], ap2->dimensions[1],
			1.0, (double *)ap2->data, lda,
			(double *)ap1->data, 1, 0.0, (double *)ret->data, 1);
	}
	else if (typenum == PyArray_CDOUBLE) {
	    cblas_zgemv(CblasRowMajor,
			CblasNoTrans,  ap2->dimensions[0], ap2->dimensions[1],
			oneD, (double *)ap2->data, lda,
			(double *)ap1->data, 1, zeroD, (double *)ret->data, 1);
	}
	else if (typenum == PyArray_FLOAT) {
	    cblas_sgemv(CblasRowMajor,
			CblasNoTrans,  ap2->dimensions[0], ap2->dimensions[1],
			1.0, (float *)ap2->data, lda,
			(float *)ap1->data, 1, 0.0, (float *)ret->data, 1);
	}
	else if (typenum == PyArray_CFLOAT) {
	    cblas_cgemv(CblasRowMajor,
			CblasNoTrans,  ap2->dimensions[0], ap2->dimensions[1],
			oneF, (float *)ap2->data, lda,
			(float *)ap1->data, 1, zeroF, (float *)ret->data, 1);
	}
    }
    else { /* (ap1->nd == 2 && ap2->nd == 2) */
	/* Matrix matrix multiplication -- Level 3 BLAS */
	lda = (ap1->dimensions[1] > 1 ? ap1->dimensions[1] : 1);
	ldb = (ap2->dimensions[1] > 1 ? ap2->dimensions[1] : 1);
	ldc = (ret->dimensions[1] > 1 ? ret->dimensions[1] : 1);
	if (typenum == PyArray_DOUBLE) {
	    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
			ap1->dimensions[0], ap2->dimensions[0], ap1->dimensions[1],
			1.0, (double *)ap1->data, lda,
			(double *)ap2->data, ldb,
			0.0, (double *)ret->data, ldc);
	}
	else if (typenum == PyArray_FLOAT) {
	    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
			ap1->dimensions[0], ap2->dimensions[0], ap1->dimensions[1],
			1.0, (float *)ap1->data, lda,
			(float *)ap2->data, ldb,
			0.0, (float *)ret->data, ldc);
	}
	else if (typenum == PyArray_CDOUBLE) {
	    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
			ap1->dimensions[0], ap2->dimensions[0], ap1->dimensions[1],
			oneD, (double *)ap1->data, lda,
			(double *)ap2->data, ldb,
			zeroD, (double *)ret->data, ldc);
	}
	else if (typenum == PyArray_CFLOAT) {
	    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
			ap1->dimensions[0], ap2->dimensions[0], ap1->dimensions[1],
			oneF, (float *)ap1->data, lda,
			(float *)ap2->data, ldb,
			zeroF, (float *)ret->data, ldc);
	}
    }
    NPY_END_ALLOW_THREADS
    Py_DECREF(ap1);
    Py_DECREF(ap2);
    return PyArray_Return(ret);

 fail:
    Py_XDECREF(ap1);
    Py_XDECREF(ap2);
    Py_XDECREF(ret);
    return NULL;
}


static char doc_vdot[] = "vdot(a,b)\nReturns the dot product of a and b for scalars and vectors\nof floating point and complex types.  The first argument, a, is conjugated.";


static PyObject *dotblas_vdot(PyObject *dummy, PyObject *args) {
    PyObject *op1, *op2;
    PyArrayObject *ap1=NULL, *ap2=NULL, *ret=NULL;
    int l;
    int typenum;
    intp dimensions[MAX_DIMS];
    PyArray_Descr *type;

    if (!PyArg_ParseTuple(args, "OO", &op1, &op2)) return NULL;

    /*
     * Conjugating dot product using the BLAS for vectors.
     * Multiplies op1 and op2, each of which must be vector.
     */

    typenum = PyArray_ObjectType(op1, 0);
    typenum = PyArray_ObjectType(op2, typenum);

    type = PyArray_DescrFromType(typenum);
    Py_INCREF(type);
    ap1 = (PyArrayObject *)PyArray_FromAny(op1, type, 0, 0, 0, NULL);
    if (ap1==NULL) {Py_DECREF(type); goto fail;}
    op1 = PyArray_Flatten(ap1, 0);
    if (op1==NULL) {Py_DECREF(type); goto fail;}
    Py_DECREF(ap1);
    ap1 = (PyArrayObject *)op1;

    ap2 = (PyArrayObject *)PyArray_FromAny(op2, type, 0, 0, 0, NULL);
    if (ap2==NULL) goto fail;
    op2 = PyArray_Flatten(ap2, 0);
    if (op2 == NULL) goto fail;
    Py_DECREF(ap2);
    ap2 = (PyArrayObject *)op2;

    if (typenum != PyArray_FLOAT && typenum != PyArray_DOUBLE &&
	typenum != PyArray_CFLOAT && typenum != PyArray_CDOUBLE) {
	if (!altered) {
	    /* need to alter dot product */
	    PyObject *tmp1, *tmp2;
	    tmp1 = PyTuple_New(0);
	    tmp2 = dotblas_alterdot(NULL, tmp1);
	    Py_DECREF(tmp1);
	    Py_DECREF(tmp2);
	}
	if (PyTypeNum_ISCOMPLEX(typenum)) {
            op1 = PyArray_Conjugate(ap1, NULL);
	    if (op1==NULL) goto fail;
	    Py_DECREF(ap1);
	    ap1 = (PyArrayObject *)op1;
	}
	ret = (PyArrayObject *)PyArray_InnerProduct((PyObject *)ap1,
						    (PyObject *)ap2);
	Py_DECREF(ap1);
	Py_DECREF(ap2);
	return PyArray_Return(ret);
    }

    if (ap2->dimensions[0] != ap1->dimensions[ap1->nd-1]) {
	PyErr_SetString(PyExc_ValueError, "vectors have different lengths");
	goto fail;
    }
    l = ap1->dimensions[ap1->nd-1];

    ret = (PyArrayObject *)PyArray_SimpleNew(0, dimensions, typenum);
    if (ret == NULL) goto fail;

    NPY_BEGIN_ALLOW_THREADS

    /* Dot product between two vectors -- Level 1 BLAS */
    if (typenum == PyArray_DOUBLE) {
	*((double *)ret->data) = cblas_ddot(l, (double *)ap1->data, 1,
					    (double *)ap2->data, 1);
    }
    else if (typenum == PyArray_FLOAT) {
	*((float *)ret->data) = cblas_sdot(l, (float *)ap1->data, 1,
					   (float *)ap2->data, 1);
    }
    else if (typenum == PyArray_CDOUBLE) {
	cblas_zdotc_sub(l, (double *)ap1->data, 1,
			(double *)ap2->data, 1, (double *)ret->data);
    }
    else if (typenum == PyArray_CFLOAT) {
	cblas_cdotc_sub(l, (float *)ap1->data, 1,
			(float *)ap2->data, 1, (float *)ret->data);
    }

    NPY_END_ALLOW_THREADS

    Py_DECREF(ap1);
    Py_DECREF(ap2);
    return PyArray_Return(ret);

 fail:
    Py_XDECREF(ap1);
    Py_XDECREF(ap2);
    Py_XDECREF(ret);
    return NULL;
}

static struct PyMethodDef dotblas_module_methods[] = {
    {"dot",  (PyCFunction)dotblas_matrixproduct, 1, doc_matrixproduct},
    {"inner",   (PyCFunction)dotblas_innerproduct,  1, doc_innerproduct},
    {"vdot", (PyCFunction)dotblas_vdot, 1, doc_vdot},
    {"alterdot", (PyCFunction)dotblas_alterdot, 1, doc_alterdot},
    {"restoredot", (PyCFunction)dotblas_restoredot, 1, doc_restoredot},
    {NULL,		NULL, 0}		/* sentinel */
};

/* Initialization function for the module */
PyMODINIT_FUNC init_dotblas(void) {
    int i;
    PyObject *d, *s;

    /* Create the module and add the functions */
    Py_InitModule3("_dotblas", dotblas_module_methods, module_doc);

    /* Import the array object */
    import_array();

    /* Initialise the array of dot functions */
    for (i = 0; i < PyArray_NTYPES; i++)
	oldFunctions[i] = NULL;

    /* alterdot at load */
    d = PyTuple_New(0);
    s = dotblas_alterdot(NULL, d);
    Py_DECREF(d);
    Py_DECREF(s);

}
