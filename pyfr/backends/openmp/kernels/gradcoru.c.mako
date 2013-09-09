# -*- coding: utf-8 -*-

<%namespace name='util' module='pyfr.backends.openmp.makoutil' />
<%include file='common.h.mako' />

static NOINLINE void
gradcoru_aux(int neles,
             ${util.arr_args('jm', [ndims, ndims], const=True)},
             ${util.arr_args('tgrad_u', [ndims, nvars])})
{
    ${util.arr_align('jm', [ndims, ndims])};
    ${util.arr_align('tgrad_u', [ndims, nvars])};

    for (int eidx = 0; eidx < neles; eidx++)
    {
        // Dereference the (transformed) gradient
    % for i, j in util.ndrange(ndims, nvars):
        ${dtype} gu${i}${j} = tgrad_u${i}${j}[eidx];
    % endfor

        // Untransform and store
    % for i, j in util.ndrange(ndims, nvars):
        tgrad_u${i}${j}[eidx] = ${' + '.join('jm{0}{2}[eidx]*gu{2}{1}'
                                             .format(i, j, k)
                                             for k in range(ndims))};
    % endfor
    }
}

void
gradcoru(int nfpts, int neles,
         const ${dtype} *jmats, ${dtype} *tgrad_u,
         int lsdj, int lsdg)
{
    #pragma omp parallel for
    for (int fidx = 0; fidx < nfpts; fidx++)
    {
        gradcoru_aux(neles,
                     ${', '.join('jmats + (fidx*{} + {})*lsdj'
                                 .format(ndims**2, i)
                                 for i in range(ndims**2))},
                     ${', '.join('tgrad_u + (({}*nfpts + fidx)*{} + {})*lsdg'
                                 .format(i, nvars, j)
                                 for i, j in util.ndrange(ndims, nvars))});
    }
}
