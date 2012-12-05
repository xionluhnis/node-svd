#include "node.h"
#include "v8.h"

extern "C" {
#include "svdlib.h"
}

using namespace v8;
using namespace node;

/**
 * Svd computation using svdlibc
 * 
 * @param args [A, dim, untranspose, debug]
 * @return the computation result
 */
Handle<Value> Svd(const Arguments& args) {
    HandleScope scope;

    // parsing the arguments
    int rows = 0, cols = 0;
    int dim = 0; // number of dimensions, 0=all by default
    bool untranspose = false; // whether to untranspose the result
    Local<Array> m = Array::New(0);
    bool merr = false;
    switch (args.Length()) {
        case 4: // with debug level
            if (!args[3]->IsNumber()) {
                return ThrowException(Exception::Error(String::New("Debug type must be a number!")));
            }
            SVDVerbosity = long(args[3]->NumberValue());
        case 3: // with untranspose
            if (!args[2]->IsBoolean()) {
                return ThrowException(Exception::Error(String::New("Untranspose must be a boolean!")));
            }
            untranspose = args[2]->BooleanValue();
        case 2: // with dimension
            if (!args[1]->IsNumber()) {
                return ThrowException(Exception::Error(String::New("Dimension must be a number!")));
            }
            dim = int(args[1]->NumberValue());
        case 1: // only A
            if (!args[0]->IsArray()) {
                merr = true;
            } else {
                m = args[0].As<Array > ();
                rows = m->Length();
                if (rows == 0) {
                    merr = true;
                } else {
                    Local<Value> v0 = m->Get(Number::New(0));
                    if (!v0->IsArray()) merr = true;
                    else {
                        Local<Array> r0 = v0.As<Array > ();
                        cols = r0->Length();
                        if (cols == 0) merr = true;
                    }
                }
            }
            if (merr) {
                return ThrowException(Exception::Error(String::New("First argument must be a non-empty matrix of number!")));
            }
            break;
        default:
            return ThrowException(Exception::Error(String::New("Requires at least a matrix argument!")));
    }

    // 1 = create matrix and populate it
    DMat D = svdNewDMat(rows, cols);
    for (int y = 0; y < rows; ++y) {
        Local<Array> row = m->Get(y).As<Array > ();
        for (int x = 0; x < cols; ++x) {
            D->value[y][x] = row->Get(x)->NumberValue();
        }
    }
    // sparse matrix generation
    SMat A = svdConvertDtoS(D);
    svdFreeDMat(D); // we can release the dense matrix

    // 2 = svd !
    SVDRec rec = svdLAS2A(A, dim);
    svdFreeSMat(A);

    // 3 = wrap result
    Local<Object> res = Object::New();
    // Dimension
    res->Set(String::New("d"), Number::New(dim = rec->d));
    // Ut vector
    Local<Array> U;
    if (untranspose) {
        U = Array::New(rows);
        for (int y = 0; y < rows; ++y) {
            Local<Array> row = Array::New(dim);
            for (int x = 0; x < dim; ++x) {
                row->Set(x, Number::New(rec->Ut->value[x][y]));
            }
            U->Set(y, row);
        }
    } else {
        U = Array::New(dim);
        for (int x = 0; x < dim; ++x) {
            Local<Array> trow = Array::New(rows);
            for (int y = 0; y < rows; ++y) {
                trow->Set(y, Number::New(rec->Ut->value[x][y]));
            }
            U->Set(x, trow);
        }
    }
    res->Set(String::New("U"), U);
    // Singular values
    Local<Array> S = Array::New(rec->d);
    for (int s = 0; s < rec->d; ++s) S->Set(s, Number::New(rec->S[s]));
    res->Set(String::New("S"), S);
    // Vt vector
    Local<Array> V;
    if (untranspose) {
        V = Array::New(cols);
        for (int y = 0; y < cols; ++y) {
            Local<Array> row = Array::New(dim);
            for (int x = 0; x < dim; ++x) {
                row->Set(x, Number::New(rec->Vt->value[x][y]));
            }
            V->Set(y, row);
        }
    } else {
        V = Array::New(dim);
        for (int x = 0; x < dim; ++x) {
            Local<Array> trow = Array::New(cols);
            for (int y = 0; y < cols; ++y) {
                trow->Set(y, Number::New(rec->Vt->value[x][y]));
            }
            V->Set(x, trow);
        }
    }
    res->Set(String::New("V"), V);
    // release result memory finally
    svdFreeSVDRec(rec);

    return scope.Close(res);
}

void Init(Handle<Object> target) {
    target->Set(String::NewSymbol("svd"), FunctionTemplate::New(Svd)->GetFunction());
}

NODE_MODULE(svd, Init)