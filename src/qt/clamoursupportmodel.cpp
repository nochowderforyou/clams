#include "clamoursupportmodel.h"

#include "guiutil.h"

ClamourSupportModel::ClamourSupportModel(QWidget *parent) :
    QAbstractTableModel(parent)
{
    columns << tr("Petition ID") << tr("Support");
}

int ClamourSupportModel::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return support.size();
}

int ClamourSupportModel::columnCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent);
    return columns.length();
}

QVariant ClamourSupportModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    std::pair<std::string, int> petition = support.at(index.row());
    if (role == Qt::DisplayRole || role == Qt::EditRole)
    {
        switch (index.column())
        {
        case 0:
            return QString::fromStdString(petition.first);
        case 1:
            return petition.second;
        }
    }
    return QVariant();
}

QVariant ClamourSupportModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (orientation != Qt::Horizontal)
        return QVariant();
    if (role == Qt::DisplayRole)
    {
        return columns.at(section);
    }
    return QVariant();
}

Qt::ItemFlags ClamourSupportModel::flags(const QModelIndex &index) const
{
    if (!index.isValid())
        return 0;
    return Qt::ItemFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
}

void ClamourSupportModel::clear()
{
    beginResetModel();
    support.clear();
    endResetModel();
}

void ClamourSupportModel::setSupport(std::vector<std::pair<std::string, int> > newSupport)
{
    beginResetModel();
    support = newSupport;
    endResetModel();
}
