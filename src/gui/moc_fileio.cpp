/****************************************************************************
** Meta object code from reading C++ file 'fileio.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "fileio.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'fileio.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_FileIO_t {
    QByteArrayData data[34];
    char stringdata0[351];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_FileIO_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_FileIO_t qt_meta_stringdata_FileIO = {
    {
QT_MOC_LITERAL(0, 0, 6), // "FileIO"
QT_MOC_LITERAL(1, 7, 13), // "sourceChanged"
QT_MOC_LITERAL(2, 21, 0), // ""
QT_MOC_LITERAL(3, 22, 6), // "source"
QT_MOC_LITERAL(4, 29, 15), // "grSourceChanged"
QT_MOC_LITERAL(5, 45, 8), // "grSource"
QT_MOC_LITERAL(6, 54, 18), // "inputSourceChanged"
QT_MOC_LITERAL(7, 73, 11), // "inputSource"
QT_MOC_LITERAL(8, 85, 16), // "linecountChanged"
QT_MOC_LITERAL(9, 102, 4), // "int&"
QT_MOC_LITERAL(10, 107, 9), // "linecount"
QT_MOC_LITERAL(11, 117, 14), // "homeDirChanged"
QT_MOC_LITERAL(12, 132, 7), // "homeDir"
QT_MOC_LITERAL(13, 140, 17), // "workingDirChanged"
QT_MOC_LITERAL(14, 158, 10), // "workingDir"
QT_MOC_LITERAL(15, 169, 11), // "modeChanged"
QT_MOC_LITERAL(16, 181, 4), // "mode"
QT_MOC_LITERAL(17, 186, 5), // "error"
QT_MOC_LITERAL(18, 192, 3), // "msg"
QT_MOC_LITERAL(19, 196, 9), // "setSource"
QT_MOC_LITERAL(20, 206, 12), // "setLinecount"
QT_MOC_LITERAL(21, 219, 10), // "setHomeDir"
QT_MOC_LITERAL(22, 230, 11), // "setgrSource"
QT_MOC_LITERAL(23, 242, 14), // "setinputSource"
QT_MOC_LITERAL(24, 257, 13), // "setWorkingDir"
QT_MOC_LITERAL(25, 271, 7), // "setMode"
QT_MOC_LITERAL(26, 279, 4), // "read"
QT_MOC_LITERAL(27, 284, 4), // "data"
QT_MOC_LITERAL(28, 289, 7), // "read_gr"
QT_MOC_LITERAL(29, 297, 10), // "read_input"
QT_MOC_LITERAL(30, 308, 15), // "startSimulation"
QT_MOC_LITERAL(31, 324, 5), // "data2"
QT_MOC_LITERAL(32, 330, 14), // "killSimulation"
QT_MOC_LITERAL(33, 345, 5) // "write"

    },
    "FileIO\0sourceChanged\0\0source\0"
    "grSourceChanged\0grSource\0inputSourceChanged\0"
    "inputSource\0linecountChanged\0int&\0"
    "linecount\0homeDirChanged\0homeDir\0"
    "workingDirChanged\0workingDir\0modeChanged\0"
    "mode\0error\0msg\0setSource\0setLinecount\0"
    "setHomeDir\0setgrSource\0setinputSource\0"
    "setWorkingDir\0setMode\0read\0data\0read_gr\0"
    "read_input\0startSimulation\0data2\0"
    "killSimulation\0write"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_FileIO[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      21,   14, // methods
       7,  184, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       8,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,  119,    2, 0x06 /* Public */,
       4,    1,  122,    2, 0x06 /* Public */,
       6,    1,  125,    2, 0x06 /* Public */,
       8,    1,  128,    2, 0x06 /* Public */,
      11,    1,  131,    2, 0x06 /* Public */,
      13,    1,  134,    2, 0x06 /* Public */,
      15,    1,  137,    2, 0x06 /* Public */,
      17,    1,  140,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      19,    1,  143,    2, 0x0a /* Public */,
      20,    1,  146,    2, 0x0a /* Public */,
      21,    1,  149,    2, 0x0a /* Public */,
      22,    1,  152,    2, 0x0a /* Public */,
      23,    1,  155,    2, 0x0a /* Public */,
      24,    1,  158,    2, 0x0a /* Public */,
      25,    1,  161,    2, 0x0a /* Public */,

 // methods: name, argc, parameters, tag, flags
      26,    1,  164,    2, 0x02 /* Public */,
      28,    1,  167,    2, 0x02 /* Public */,
      29,    1,  170,    2, 0x02 /* Public */,
      30,    2,  173,    2, 0x02 /* Public */,
      32,    1,  178,    2, 0x02 /* Public */,
      33,    1,  181,    2, 0x02 /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, QMetaType::QString,    5,
    QMetaType::Void, QMetaType::QString,    7,
    QMetaType::Void, 0x80000000 | 9,   10,
    QMetaType::Void, QMetaType::QString,   12,
    QMetaType::Void, QMetaType::QString,   14,
    QMetaType::Void, QMetaType::QString,   16,
    QMetaType::Void, QMetaType::QString,   18,

 // slots: parameters
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, 0x80000000 | 9,   10,
    QMetaType::Void, QMetaType::QString,   12,
    QMetaType::Void, QMetaType::QString,    5,
    QMetaType::Void, QMetaType::QString,    7,
    QMetaType::Void, QMetaType::QString,   14,
    QMetaType::Void, QMetaType::QString,   16,

 // methods: parameters
    QMetaType::QString, QMetaType::QString,   27,
    QMetaType::QString, QMetaType::QString,   27,
    QMetaType::QString, QMetaType::QString,   27,
    QMetaType::Int, QMetaType::QString, QMetaType::QString,   27,   31,
    QMetaType::Int, QMetaType::QString,   27,
    QMetaType::Bool, QMetaType::QString,   27,

 // properties: name, type, flags
       3, QMetaType::QString, 0x00495103,
       5, QMetaType::QString, 0x00495003,
       7, QMetaType::QString, 0x00495003,
      10, QMetaType::Int, 0x00495103,
      12, QMetaType::QString, 0x00495103,
      14, QMetaType::QString, 0x00495103,
      16, QMetaType::QString, 0x00495103,

 // properties: notify_signal_id
       0,
       1,
       2,
       3,
       4,
       5,
       6,

       0        // eod
};

void FileIO::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FileIO *_t = static_cast<FileIO *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->sourceChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: _t->grSourceChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 2: _t->inputSourceChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 3: _t->linecountChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->homeDirChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: _t->workingDirChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: _t->modeChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: _t->error((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 8: _t->setSource((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 9: _t->setLinecount((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: _t->setHomeDir((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 11: _t->setgrSource((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 12: _t->setinputSource((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 13: _t->setWorkingDir((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 14: _t->setMode((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 15: { QString _r = _t->read((*reinterpret_cast< const QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< QString*>(_a[0]) = std::move(_r); }  break;
        case 16: { QString _r = _t->read_gr((*reinterpret_cast< const QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< QString*>(_a[0]) = std::move(_r); }  break;
        case 17: { QString _r = _t->read_input((*reinterpret_cast< const QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< QString*>(_a[0]) = std::move(_r); }  break;
        case 18: { int _r = _t->startSimulation((*reinterpret_cast< const QString(*)>(_a[1])),(*reinterpret_cast< const QString(*)>(_a[2])));
            if (_a[0]) *reinterpret_cast< int*>(_a[0]) = std::move(_r); }  break;
        case 19: { int _r = _t->killSimulation((*reinterpret_cast< const QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< int*>(_a[0]) = std::move(_r); }  break;
        case 20: { bool _r = _t->write((*reinterpret_cast< const QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = std::move(_r); }  break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::sourceChanged)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::grSourceChanged)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::inputSourceChanged)) {
                *result = 2;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(int & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::linecountChanged)) {
                *result = 3;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::homeDirChanged)) {
                *result = 4;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::workingDirChanged)) {
                *result = 5;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::modeChanged)) {
                *result = 6;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::error)) {
                *result = 7;
                return;
            }
        }
    }
#ifndef QT_NO_PROPERTIES
    else if (_c == QMetaObject::ReadProperty) {
        FileIO *_t = static_cast<FileIO *>(_o);
        Q_UNUSED(_t)
        void *_v = _a[0];
        switch (_id) {
        case 0: *reinterpret_cast< QString*>(_v) = _t->source(); break;
        case 1: *reinterpret_cast< QString*>(_v) = _t->grSource(); break;
        case 2: *reinterpret_cast< QString*>(_v) = _t->inputSource(); break;
        case 3: *reinterpret_cast< int*>(_v) = _t->linecount(); break;
        case 4: *reinterpret_cast< QString*>(_v) = _t->homeDir(); break;
        case 5: *reinterpret_cast< QString*>(_v) = _t->workingDir(); break;
        case 6: *reinterpret_cast< QString*>(_v) = _t->mode(); break;
        default: break;
        }
    } else if (_c == QMetaObject::WriteProperty) {
        FileIO *_t = static_cast<FileIO *>(_o);
        Q_UNUSED(_t)
        void *_v = _a[0];
        switch (_id) {
        case 0: _t->setSource(*reinterpret_cast< QString*>(_v)); break;
        case 1: _t->setgrSource(*reinterpret_cast< QString*>(_v)); break;
        case 2: _t->setinputSource(*reinterpret_cast< QString*>(_v)); break;
        case 3: _t->setLinecount(*reinterpret_cast< int*>(_v)); break;
        case 4: _t->setHomeDir(*reinterpret_cast< QString*>(_v)); break;
        case 5: _t->setWorkingDir(*reinterpret_cast< QString*>(_v)); break;
        case 6: _t->setMode(*reinterpret_cast< QString*>(_v)); break;
        default: break;
        }
    } else if (_c == QMetaObject::ResetProperty) {
    }
#endif // QT_NO_PROPERTIES
}

const QMetaObject FileIO::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_FileIO.data,
      qt_meta_data_FileIO,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *FileIO::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *FileIO::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_FileIO.stringdata0))
        return static_cast<void*>(const_cast< FileIO*>(this));
    return QObject::qt_metacast(_clname);
}

int FileIO::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 21)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 21;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 21)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 21;
    }
#ifndef QT_NO_PROPERTIES
   else if (_c == QMetaObject::ReadProperty || _c == QMetaObject::WriteProperty
            || _c == QMetaObject::ResetProperty || _c == QMetaObject::RegisterPropertyMetaType) {
        qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyDesignable) {
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyScriptable) {
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyStored) {
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyEditable) {
        _id -= 7;
    } else if (_c == QMetaObject::QueryPropertyUser) {
        _id -= 7;
    }
#endif // QT_NO_PROPERTIES
    return _id;
}

// SIGNAL 0
void FileIO::sourceChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void FileIO::grSourceChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void FileIO::inputSourceChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void FileIO::linecountChanged(int & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void FileIO::homeDirChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void FileIO::workingDirChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void FileIO::modeChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void FileIO::error(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
