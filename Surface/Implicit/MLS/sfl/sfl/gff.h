/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   gff.h
 * Created: Thu Apr 11 23:29:55 2002
 * Authors: Tim Weyrich <weyrich@inf.ethz.ch>
 *
 * Copyright (c) 2001--2003, Computer Graphics Lab, ETH Zurich
 *
 * This file is part of the Pointshop3D system.
 * See http://www.pointshop3d.com/ for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307 USA
 *
 * Contact info@pointshop3d.com if any conditions of this
 * licensing are not clear to you.
 *
 * ---------------------------------------------------------------
 *
 * $Log: gff.h,v $
 * Revision 1.7  2003/11/19 13:04:01  weyrich
 * added Pointshop3D headers
 *
 * Revision 1.6  2003/11/17 15:05:52  weyrich
 * new API allowing for more than 32 surfel properties
 *
 * Revision 1.5  2003/03/20 14:52:51  weyrich
 * smaller fixes of implicit coordinate system transformations
 *
 * Revision 1.4  2002/11/18 13:09:12  weyrich
 * added markers and tangential vectors
 *
 * Revision 1.3  2002/05/17 13:58:16  weyrich
 * PointShop3D integration nearly finished
 *
 * Revision 1.2  2002/05/08 09:20:36  weyrich
 * added PI_PACKED/pragma pack to disc structures
 *
 * Revision 1.1.1.1  2002/05/08 08:58:47  weyrich
 * surfel file format
 */

#ifndef __GFF_H__
#define __GFF_H__

#include <pitypes.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

  /*
   *  Win32 Stuff
   */

#ifdef _WIN32
  /* The following ifdef block is the standard way of creating macros which make exporting
   * from a DLL simpler. All files within this DLL are compiled with the GFF_EXPORTS
   * symbol defined on the command line. this symbol should not be defined on any project
   * that uses this DLL. This way any other project whose source files include this file see
   * GFF_API functions as being imported from a DLL, wheras this DLL sees symbols
   * defined with this macro as being exported.
   */
#  ifdef GFF_EXPORTS
#    define GFF_API __declspec(dllexport)
#    include "../stdafx.h"
#  else
#    define GFF_API __declspec(dllimport)
#  endif
#else
#  define GFF_API
#endif

  /*
   *  File Format Definitions
   */

#define  GFF_MAGIC    "GFLF"
#define  GFF_VERSION  0x0000

  /* GFFTAG_LINK is optional and overrides the guess, that the next
     directory starts after the GFFTAG_DATA data.
  */

#define  GFFTAG_UPLINK        0x0000
#define  GFFTAG_LINK          0x0001
#define  GFFTAG_DOWNLINK      0x0002
#define  GFFTAG_NUM_CHILDREN  0x0003

#define  GFF_INTERNAL_MAXTAG  0x0004

  extern const char  *gff_tag_str[GFF_INTERNAL_MAXTAG];

#define  GFF_P3D_STARTTAG     0x0100
#define  GFF_P3D_MAXTAG       0x0500

#define  GFFTYPE_UINT8        0
#define  GFFTYPE_INT8         1
#define  GFFTYPE_UINT16       2
#define  GFFTYPE_INT16        3
#define  GFFTYPE_UINT32       4
#define  GFFTYPE_INT32        5
#define  GFFTYPE_FLOAT        6
#define  GFFTYPE_DOUBLE       7
#define  GFFTYPE_ASCII        8
#define  GFFTYPE_APPLICATION  9 /*!< \brief For convenience. */

  /* To assign LDB data: */
#define  GFFTYPE_LDB      10    /*!< Can only be produced by gff_start_data(). */

  /* For internal use only: */
#define  GFFTYPE_OFFSET   11
#define  GFFTYPE_2OFFSET  12

#define  GFFTYPE_NUM_TYPES  13

#define  GFF_UNKNOWN_SIZE  0xffffffffUL

  const uint32  gff_type_size[GFFTYPE_NUM_TYPES] =
    {
      1, 1, 2, 2, 4, 4, 4, 8, 1, 1, 1, 4, 8
    };

  extern const char  *gff_type_str[GFFTYPE_NUM_TYPES];

  typedef void gff_directory;

#ifdef _WIN32
#  pragma pack(push, __pre_gff_header__)
#  pragma pack(4)
#endif

  typedef struct
  {
    uint32  gff_magic;          /*!< Contains #GFF_MAGIC. */
    uint32  gff_version;        /*!< GFF version: 16bit major, 16bit minor version. */
    uint32  inst_magic;         /*!< Format instance dependent magic. */
    uint32  inst_version;       /*!< Format instance version: 16bit major, 16bit minor version. */
    uint32  root;               /*!< Offset of root node. */
  }
  gff_header PI_PACKED;

#ifdef _WIN32
#  pragma pack(pop, __pre_gff_header__)
#endif

#define  GFF_HEADER_SIZE  20    /*!< \brief Used for padding check. */

  typedef struct
  {
    int              fd;
    struct gff_dir  *root_dir;
    struct gff_dir  *active_dir;
    int              seekable_p;
    int              read_p;

    gff_header       hd;

    struct gff_dir  *currently_read_dir;
    uint32           next_dir_to_read_offs;
    struct gff_dir  *next_dir_to_read_parent;
  }
  gff_file;

  /*! \brief Open a GFF file. */
  GFF_API gff_file  *gff_open(const char *fname, const char *mode, const char magic[4], uint32 version);
  GFF_API int        gff_close(gff_file *file);

  /*! \brief Open a GFF file from a file descriptor. */
  GFF_API gff_file  *gff_fdopen(int fd, const char *mode, const char magic[4], uint32 version);

  /*! \brief Low-level, signal-safe block read. */
  GFF_API int32  gff_read(gff_file *file, void *buf, uint32 size);
  /*! \brief Low-level, signal-safe block write. */
  GFF_API int32  gff_write(gff_file *file, const void *buf, uint32 size);

  /*
   *  Read mode
   */

  GFF_API uint32  gff_get_magic(const gff_file *file);
  GFF_API uint32  gff_get_version(const gff_file *file);
  GFF_API int     gff_is_seekable(const gff_file *file);
  GFF_API void    gff_assume_not_seekable(gff_file *file);

  GFF_API gff_directory  *gff_read_root_dir(gff_file *file);
  GFF_API gff_directory  *gff_read_next_dir(gff_file *file, gff_directory *cur_dir);
  GFF_API gff_directory  *gff_read_child_dir(gff_file *file, gff_directory *cur_dir, int child_idx);

  GFF_API uint32                gff_dir_get_num_children(const gff_directory *dir);
  GFF_API uint32                gff_dir_get_num_tags(const gff_directory *dir);
  GFF_API int                   gff_dir_get_tag_list(const gff_directory *dir, int *tag_ids);
  GFF_API const gff_directory  *gff_dir_get_parent_dir(const gff_directory *dir);

  GFF_API int  gff_query_field(const gff_directory *dir, int tag_id, int *field_type, int *num_elements);
  GFF_API int  gff_get_field(const gff_directory *dir, int tag_id, int field_type, int expected_size, ...);
  GFF_API int  gff_v_get_field(const gff_directory *dir, int tag_id, int field_type, int expected_size, va_list ap);
  GFF_API int  gff_get_field_arr(const gff_directory *dir, int tag_id, int field_type, int expected_size, void *dest);
  GFF_API const void  *gff_get_field_arr_p(const gff_directory *dir, int tag_id, int field_type);

  GFF_API int  gff_seek_field_data(const gff_file *file, const gff_directory *dir, int tag_id);

  /*
   *  Write mode
   */

  GFF_API int  gff_write_field(gff_file *file, int tag_id, int field_type, int num_elements, ...);
  GFF_API int  gff_v_write_field(gff_file *file, int tag_id, int field_type, int num_elements, va_list ap);
  GFF_API int  gff_write_field_arr(gff_file *file, int tag_id, int field_type, int num_elements, const void *data);
  GFF_API int  gff_write_field_cstring(gff_file *file, int tag_id, const char *str);

  GFF_API int  gff_write_downlink(gff_file *file, uint32 num_children);
  GFF_API int  gff_write_uplink(gff_file *file);

  GFF_API int  gff_write_data_field(gff_file *file, int tag_id, int num_bytes);
  GFF_API int  gff_break_dir(gff_file *file);

#ifdef __cplusplus
}
#endif

#endif

/* Some Emacs-Hints -- please don't remove:
 *
 *  Local Variables:
 *  tab-width:4
 *  End:
 */
